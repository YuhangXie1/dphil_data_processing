"""
plate_reader_analysis_diya_jbl137_wm.py
=======================================

Step 2 of the plate-reader pipeline.

Loads the ``*_extracted_combined.csv`` file made by
``plate_reader_extraction.py``, attaches per-well experimental metadata from
``experiment_input.xlsx`` (cell strain, media, green / red light intensity),
computes fluorescence normalised to OD600, averages replicate wells, and plots
the results.

Pipeline overview
-----------------
    t0.xlsx ... t24.xlsx  --(plate_reader_extraction.py)-->  *_extracted_combined.csv
                                                                     |
    experiment_input.xlsx (plate map) -----------------------------+
                                                                     v
                                                     this script --> <folder>/figures/*.png

Inputs
------
``*_extracted_combined.csv``
    Two-row column header ``(measurement, timestamp)``, one row per well.
``experiment_input.xlsx``
    Plate map in the same layout as the downstairs plate reader export: a
    ``"Label: <field>"`` row followed by 8 rows (A-H) x 12 columns of values.
    Expected labels: ``cells``, ``media``, ``green_intensity``,
    ``red_intensity``.

Key data structures
-------------------
``sorted_data_df`` ("raw", long / tidy format) - one row per
well x measurement x timepoint, columns::

    well, cells, media, green_intensity, red_intensity,
    timestamp, timepoint, time, measurement, value

Measurements include the raw channels (``OD600``, ``GFP 395nm``,
``GFP 488nm``) and derived ratios (``GFP 395nm/OD600``, ``GFP 488nm/OD600``).

``summary_df`` - replicate wells with identical metadata averaged together,
one row per condition x measurement x timepoint, columns::

    cells, media, green_intensity, red_intensity, timestamp, time, timepoint,
    measurement, min_timestamp, mean, std, count

Plot styling
------------
Each plot encodes three variables:
    * marker shape / marker colour  -> cell type    (markerstyle_map, markercolor_map)
    * line style                    -> media        (linestyle_map)
    * line colour                   -> green light  (line_color_map)
A value missing from a map falls back to the ``default_*`` style in
:class:`PlotConfig` (see :func:`style_for`). Use ``plot_exclude`` to hide
conditions. Legends are built by :func:`add_legends` and list only the
conditions that are actually plotted, stacked so they never overlap.

Tip: on each run :func:`print_style_templates` prints copy-paste-ready versions
of all these dictionaries based on the plate map.

Usage
-----
1. Run ``plate_reader_extraction.py`` first.
2. Edit the USER VARIABLES block and the file paths / plot calls in
   :func:`main`.
3. Run ``python plate_reader_analysis_diya_jbl137_wm.py``.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re
import warnings
from scipy import stats as sps
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from enum import Enum  # TIDY: unused, can be removed
from dataclasses import dataclass
from typing import Dict, Any
import seaborn as sns  # not used yet - see plotting suggestions


### USER VARIABLES ###

### HUMAN-READABLE CELL LABELS ###
# Maps the short code used in experiment_input.xlsx to a readable name.
# Currently only used for bar-chart legends in plot_comparison_bars().
# Codes ending "_<n>" are individual biological replicates (colonies).
cell_label_map = {
    "CS": "Double plasmid",
    "CS_1": "Double_plasmid_1",
    "CS_2": "Double_plasmid_2",
    "CS_3": "Double_plasmid_3",
    "CS_4": "Double_plasmid_4",
    "CS_5": "Double_plasmid_5",
    "CS_6": "Double_plasmid_6",
    "CS_7": "Double_plasmid_7",
    "CS_8": "Double_plasmid_8",
    "KS": "Single plasmid",
    "KS_1": "Single_plasmid_1",
    "KS_2": "Single_plasmid_2",
    "KS_3": "Single_plasmid_3",
    "KS_4": "Single_plasmid_4",
    "KS_5": "Single_plasmid_5",
    "KS_6": "Single_plasmid_6",
    "KS_7": "Single_plasmid_7",
    "KS_8": "Single_plasmid_8",
    "WT": "Wildtype",
    "WT_1": "Wildtype_1",
    "WT_2": "Wildtype_2",
    "WT_3": "Wildtype_3",
    "media": "media",
    "none": "none"
}

# Marker shape per cell type (matplotlib marker codes: "o" circle, "s" square,
# "x" cross, "^" triangle). Cell types without an entry (e.g. WT_1) use
# PlotConfig.default_markerstyle. Legends only show plotted cell types.
markerstyle_map = {
    "CS": "o",
    #"CS_1": "o",
    #"CS_2": "o",
    #"CS_3": "o",
    #"CS_4": "o",
    #"CS_5": "o",
    #"CS_6": "o",
    #"CS_7": "o",
    #"CS_8": "o",
    "KS": "s",
    #"KS_1": "s",
    #"KS_2": "s",
    #"KS_3": "s",
    #"KS_4": "s",
    #"KS_5": "s",
    #"KS_6": "s",
    #"KS_7": "s",
    #"KS_8": "s",
    "WT": "x",
    #"WT_1": "x",
    #"WT_2": "x",
    #"WT_3": "x",
    "media": "^",
    "none": "o",
}

# Marker fill/edge colour per cell type. Missing entries use
# PlotConfig.default_markercolor.
markercolor_map = {
    "CS": "blue",
    "CS_1": "blue",
    "CS_2": "blue",
    "CS_3": "blue",
    "CS_4": "blue",
    "CS_5": "blue",
    "CS_6": "blue",
    "CS_7": "blue",
    "CS_8": "blue",
    "KS": "green",
    "KS_1": "green",
    "KS_2": "green",
    "KS_3": "green",
    "KS_4": "green",
    "KS_5": "green",
    "KS_6": "green",
    "KS_7": "green",
    "KS_8": "green",
    "WT": "black",
    "WT_1": "black",
    "WT_2": "black",
    "WT_3": "black",
    "media": "black",
    "none": "black",
}

# Line style per media type ("solid", "dashed", "dotted", "dashdot").
linestyle_map = {
    "M9": "solid",
}

# Line colour per green light intensity, as (R, G, B) tuples in 0-1.
# Here: red (no green light) -> green (max green light).
line_color_map = {
    0.0: (1.000, 0.000, 0.000),
    0.028: (0.750, 0.250, 0.000),
    0.28: (0.500, 0.500, 0.000),
    0.56: (0.250, 0.750, 0.000),
    1.4: (0.000, 1.000, 0.000),
}

# Conditions to HIDE from plots. Uncomment a value to exclude it; leave it
# commented to show it. Values must match experiment_input.xlsx exactly
# (intensities are compared as floats).
plot_exclude = {
    "cells": [
         #"CS",
         'CS_1',
         'CS_2',
         'CS_3',
         "CS_4",
         "CS_5",
         "CS_6",
         "CS_7",
         "CS_8",
         #"KS",
         'KS_1',
         'KS_2',
         'KS_3',
         "KS_4",
         "KS_5",
         "KS_6",
         "KS_7",
         "KS_8",
         "WT",
        # 'WT_1',
        # 'WT_2',
        # 'WT_3',
         'media',
        "none",
    ],
    "media": [
        # 'M9',
    ],
    "green_intensity": [
        # 0,
        # 0.028,
        # 0.28,
        # 0.56,
        # 1.4,
    ],
    "red_intensity": [
        # 2.8,
    ],
}


### MAIN ###

def main(cell_label_map, markerstyle_map, markercolor_map, linestyle_map, line_color_map, plot_exclude):
    """
    Load the data, apply the user styling, and produce the chosen plots.

    This is the "lab notebook" area of the script: edit the file paths and
    uncomment / add plot calls here for each experiment.

    :param cell_label_map: short cell code -> readable name (currently unused here)
    :param markerstyle_map: cell type -> marker shape
    :param markercolor_map: cell type -> marker colour
    :param linestyle_map: media -> line style
    :param line_color_map: green intensity -> RGB tuple
    :param plot_exclude: conditions to hide, see the USER VARIABLES block
    """

    # Combined CSV produced by plate_reader_extraction.py.
    # Figures are saved to "<this folder>/figures/".
    filepath = "26-11-19 Kirill M9 exp 2/26-11-19 Kirill M9 exp 2_diya_extracted_combined.csv"
    DefaultConfig.save_file_path = filepath.split("/")[0]

    # Plate map (which cells/media/light are in which well).
    # Prints copy-paste templates for the style dictionaries above.
    # TIDY: path is relative to the working directory, not to `filepath`'s folder.
    plate_map = load_plate_metadata("experiment_input.xlsx")
    print_style_templates(plate_map)

    # Long-format data + replicate means/SDs
    sorted_data_df, summary_df = load_data(filepath, plate_map)

    # Copy user styling into the shared config used by all plot functions.
    # Note: this overrides the automatic line_color_map built inside load_data().
    DefaultConfig.plot_exclude = plot_exclude
    DefaultConfig.markerstyle_map = markerstyle_map
    DefaultConfig.markercolor_map = markercolor_map
    DefaultConfig.linestyle_map = linestyle_map
    DefaultConfig.line_color_map = line_color_map

    # --- Timecourse plots ----------------------------------------------------
    #plot_timecourse(sorted_data_df, summary_df, "OD600", "average", title_extra= "YX comparison - no JBL001", save_image = True)
    plot_timecourse(sorted_data_df, summary_df, "OD600", "average", title_extra= "all intensities (no WT or media) ", save_image = True)
    plot_timecourse(sorted_data_df, summary_df, "GFP 395nm/OD600", "average", title_extra= "all intensities (no WT or media) ", save_image = True)


    # --- Dose-response (by intensity) plots -----------------------------------
    # TIDY: the blocks below re-define plot_exclude / line_color_map but every
    # plot call that would use them is commented out, so they currently have
    # no effect. Consider deleting or moving to a separate "recipes" section.
    #plot_by_intensity(sorted_data_df, summary_df, "OD600", "average", 5, title_extra= "t12 YX comparison", save_image = True)
    plot_exclude = {
            "cells":[
                    "none",
                    #"media",
                    #"WT",
                    #"YX001",
                    #"YX002",

                    ],
            "media":[],
            #"green_intensity":[1.4,0.56,0.028],
            #"green_intensity":[2.8,1.4,0.56,0.028],
            "green_intensity":[],
            "red_intensity":[],
        }
    DefaultConfig.plot_exclude = plot_exclude
    # Timepoint is a 0-based index of the read: 0 = first read (t0),
    # 4 = fifth read (t12), 5 = sixth read (t15), ...
    #plot_by_intensity(sorted_data_df, summary_df, "GFP 395nm/OD600", "average", 4, title_extra= "t12 - inc media and WT", save_image = True)
    #plot_by_intensity(sorted_data_df, summary_df, "GFP 395nm/OD600", "average", 5, title_extra= "t15 - inc media and WT", save_image = True)

    # --- Double vs single plasmid comparison ---------------------------------
    line_color_map = {1.4: (0.0,1.0,0.0),
                    0.56: (0.3,0.7,0.0),
                    0.28:(0.5,0.5,0.0),
                    0.028:(0.7,0.3,0.0),
                    0.0:(1.0,0.0,0.0),
        }

    plot_exclude = {
            "cells":[
                    "none",
                    "media",
                    #"WT",
                    #"YX001",
                    #"YX002",

                    ],
            "media":[],
            "green_intensity":[0.28,0.56,0.028],
            #"green_intensity":[2.8,1.4,0.56,0.028],
            #"green_intensity":[],
            "red_intensity":[],
        }
    DefaultConfig.line_color_map = line_color_map
    #plot_timecourse(sorted_data_df, summary_df, "GFP 395nm/OD600", "average", title_extra= "Double vs Single comparison", save_image = False)

    DefaultConfig.plot_exclude = plot_exclude
    #plot_timecourse(sorted_data_df, summary_df, "GFP 395nm/OD600", "average", title_extra= "Double vs Single comparison, green vs red only", save_image = False)

    # --- Bar charts with CS vs KS statistics ----------------------------------
    # Bars use their own `conditions` list, so reset exclusions to "none" only.
    DefaultConfig.plot_exclude = {"cells": ["none"], "media": [], "green_intensity": [], "red_intensity": []}
    t12 = timepoint_at_hours(sorted_data_df, 12)   # read index nearest 12 h

    # Across intensities at ~12 h and ~24 h (one panel each), raw values
    #plot_comparison_bars(sorted_data_df, "GFP 395nm/OD600", across="intensity",
    #                     timepoints=[t12, timepoint_at_hours(sorted_data_df, 24)],
    #                     title_extra="CS vs KS", save_image=True)

    # Same, as fold change vs WT
    #plot_comparison_bars(sorted_data_df, "GFP 395nm/OD600", across="intensity", timepoints=t12,
    #                     conditions=("CS", "KS", "WT"), normalise="WT", p_format="both",
    #                     title_extra="CS vs KS", save_image=True)

    # Across time at the highest and lowest green intensity
    #plot_comparison_bars(sorted_data_df, "GFP 395nm/OD600", across="timepoint", intensities=[0, 1.4],
    #                     conditions=("CS", "KS", "WT"), title_extra="CS vs KS", save_image=True)

    # Same figure with two-way ANOVA + post-hoc contrasts instead of t-tests
    #plot_comparison_bars(sorted_data_df, "GFP 395nm/OD600", across="intensity",
    #                     timepoints=[t12, timepoint_at_hours(sorted_data_df, 24)],
    #                     test="anova", title_extra="CS vs KS ANOVA", save_image=True)
    # ...ANOVA on the responsive range only (bars for 0 and 0.028 still drawn)
    #plot_comparison_bars(sorted_data_df, "GFP 395nm/OD600", across="intensity",
    #                     timepoints=[t12, timepoint_at_hours(sorted_data_df, 24)],
    #                     test="anova", stats_levels=[0.28, 0.56, 1.4],
    #                     title_extra="CS vs KS ANOVA 0.28-1.4", save_image=True)

    # Statistics table only
    #stats = compare_conditions_stats(sorted_data_df, "GFP 395nm/OD600", "CS", "KS", timepoints=t12)
    #print(stats)


### PLOT CONFIGURATION ###

@dataclass
class PlotConfig:
    """
    Styling and output settings shared by all plotting functions.

    :ivar markerstyle_map: cell type -> matplotlib marker (e.g. ``"o"``)
    :ivar markercolor_map: cell type -> marker colour
    :ivar linestyle_map: media -> line style (e.g. ``"solid"``)
    :ivar line_color_map: green intensity -> RGB(A) colour tuple
    :ivar line_alpha_map: green intensity -> line transparency (0-1);
        only used by :func:`plot_timecourse_custom` when ``alpha_used`` is True
    :ivar plot_exclude: ``{"cells": [...], "media": [...],
        "green_intensity": [...], "red_intensity": [...]}`` values to hide
    :ivar alpha_used: whether to apply ``line_alpha_map``
    :ivar save_file_path: folder in which a ``figures/`` subfolder is created
    :ivar default_*: fallback styles used by :func:`style_for` when a value
        has no entry in the relevant map
    """
    markerstyle_map: Dict[str,str]
    markercolor_map: Dict[str,str]
    linestyle_map: Dict[str,str]
    line_color_map: Dict[Any,tuple]
    line_alpha_map: Dict[Any,float]
    plot_exclude: Dict[str, list]
    alpha_used: bool
    save_file_path: str

    default_markerstyle: str = "o"
    default_markercolor: str = "blue"
    default_linestyle: str = "solid"
    default_linecolor: str = "black"
    default_linealpha: float = 1.0


def default_plot_config(save_file_path: str = ".") -> PlotConfig:
    """
    Build a :class:`PlotConfig` with starter values.

    These defaults are from an earlier experiment (JBL001 / YX001 / YX002 in
    WM / LB media) and are normally overwritten in :func:`main`.

    :param save_file_path: folder where ``figures/`` will be created
    :return: a new PlotConfig
    """
    return PlotConfig(
        markerstyle_map = {"JBL001":"s", "YX002":"*", "YX001":"o", "media":"^"},
        markercolor_map = {"JBL001":"black", "YX002":"orange", "YX001": "blue", "media":"black"},
        linestyle_map = {"WM-met+": "solid", "WM-met-": "dashed", "LB": "solid"},
        line_color_map = {},  # built from the data in load_data()
        line_alpha_map = {},
        plot_exclude = {"cells":[], "media":[], "green_intensity":[], "red_intensity":[]},
        alpha_used = False,
        save_file_path = save_file_path
    )

# Global config object shared by every plot function (used as their default
# `config` argument). main() mutates its attributes before plotting.
DefaultConfig = default_plot_config()


### FUNCTIONS ###

def mean_std_cv(data):
    """
    Mean, standard deviation and coefficient of variation across replicates.

    TIDY: legacy helper - not called anywhere in the current pipeline
    (replaced by the groupby in :func:`load_data`).

    :param data: 2D array-like, one row per replicate
        (e.g. ``[[rep1_t0, rep1_t1, ...], [rep2_t0, ...]]``)
    :return: ``[mean, std, cv]``, each computed along axis 0 (across replicates).
        ``std`` is the population SD (ddof=0); ``cv`` is in percent.
    """
    mean = np.mean(data, axis = 0)
    std_dev = np.std(data, axis = 0)
    cv = (std_dev/mean)*100
    return [mean, std_dev, cv]


def process_data_to_mean(data_dict):
    """
    Group replicate entries by name and compute mean / SD / CV for each group.

    Replicates are identified by stripping the final ``"_<suffix>"`` from each
    key, e.g. ``"YX001_WM-met+_1"`` and ``"YX001_WM-met+_2"`` are both grouped
    under ``"YX001_WM-met+"``.

    TIDY: legacy helper - not called anywhere in the current pipeline.

    :param data_dict: ``{replicate_name: data_array}``
    :return: ``{group_name: [mean, std, cv]}``
    """
    # Collect all replicate arrays under their shared group name
    data_dict_repeats = {}
    for key, value in data_dict.items():
        #if key == "media":
        #    data_dict_repeats["media"].append(value)
        # "name_media_repeat" -> ["name_media", "repeat"]
        namemet_repeat = key.rsplit("_",1)
        if namemet_repeat[0] in data_dict_repeats.keys():
            data_dict_repeats[namemet_repeat[0]].append(value)
        else:
            data_dict_repeats[namemet_repeat[0]] = [value]

    print(data_dict_repeats)
    data_dict_processed = {}
    for key, value in data_dict_repeats.items():
        # mean_std_cv_array = [[mean_od, mean_gfp, mean_gfp/od],
        #                      [std_od,  std_gfp,  std_gfp/od],
        #                      [cv_od,   cv_gfp,   cv_gfp/od]]
        mean_std_cv_array = mean_std_cv(value)

        # value = [mean, std, cv]
        data_dict_processed[key] = mean_std_cv_array
    print(data_dict_processed)
    return data_dict_processed


def build_intensity_colour_map(green_intensities, cmap_name = None) -> Dict[Any,tuple]:
    """
    Automatically assign a line colour to each green light intensity.

    Intensities are scaled linearly between the minimum and maximum value and
    mapped onto a colormap (red -> green by default).

    :param green_intensities: iterable of green intensity values (duplicates OK)
    :param cmap_name: optional matplotlib colormap name (e.g. ``"viridis"``);
        if None, a red-to-green gradient is used
    :return: ``{intensity: (r, g, b, a)}``, or ``{}`` if only one intensity
        exists (TIDY: an empty map then causes KeyError when plotting)
    """
    green_intensities_set = set(green_intensities)
    min_intensity = min(green_intensities_set)
    max_intensity = max(green_intensities_set)

    if cmap_name is None:
        cmap = LinearSegmentedColormap.from_list("red_to_green",[(1,0,0), (0,1,0)])
    else:
        cmap = plt.get_cmap(cmap_name)

    # Cannot normalise a single value (division by zero)
    if max_intensity == min_intensity:
        print("only one green intensity")
        return {}

    # Linear scaling: min -> 0.0 (red), max -> 1.0 (green)
    norm = mcolors.Normalize(vmin=min_intensity,vmax=max_intensity)

    return {
            intensity: cmap(norm(intensity))
            for intensity in green_intensities_set
        }


def create_plot_handles(config: PlotConfig):
    """
    Build legend entries (dummy ``Line2D`` handles) for the three legends.

    SUPERSEDED by :func:`add_legends`, which only shows plotted conditions.
    No longer called by any plot function; kept for reference.

    Legends are built from the style maps, not from the plotted data, so every
    key in the maps is shown - including excluded conditions.

    :param config: plot configuration
    :return: ``(cell_handles, intensity_handles, media_handles)``
        * cell_handles: one marker per key in ``markerstyle_map``
        * intensity_handles: one coloured line per key in ``line_color_map``,
          highest intensity first
        * media_handles: one black line per key in ``linestyle_map``
    """
    # Cell type legend: markers only (linestyle "None")
    cell_handles = [
        Line2D(
            [], [],
            marker = config.markerstyle_map[cell],
            linestyle = "None",
            markerfacecolor = config.markercolor_map[cell],
            markeredgecolor = config.markercolor_map[cell],
            markersize = 6,
            label = cell
        )
        for cell in config.markerstyle_map
    ]

    # Light intensity legend: sorted high -> low
    sorted_intensities = sorted(
        config.line_color_map.keys(),
        key=float,
        reverse=True
    )

    intensity_handles = [
        Line2D(
            [], [],
            color=config.line_color_map[intensity],
            linewidth=2,
            label=intensity
        )
        for intensity in sorted_intensities
    ]

    # Media legend: black lines with each media's line style
    media_handles = [
        Line2D(
            [], [],
            color = "black",
            linestyle = config.linestyle_map[media],
            linewidth = 2,
            label = media
        )
        for media in config.linestyle_map
    ]

    return cell_handles, intensity_handles, media_handles


### STYLE, LEGEND AND OUTPUT HELPERS ###
# Shared by all plot functions. They make sure that:
#   * excluded conditions are removed in ONE place (apply_exclusions),
#   * a missing style entry falls back to the PlotConfig defaults instead of
#     raising KeyError (style_for),
#   * legends list only what is actually plotted and never overlap
#     (add_legends),
#   * figures are saved with the legends included (save_or_show).

# Default columns, legend titles and order of the legends (top to bottom)
LEGEND_SPECS = {
    "cells": "Cell type",
    "media": "Media",
    "green_intensity": "Green light intensity",
}


def style_for(config, cells=None, media=None, green_intensity=None):
    """
    Look up plot styles for one condition, with fallbacks.

    Any value missing from a style map gets the matching ``default_*`` value
    from :class:`PlotConfig` (e.g. a cell type not in ``markerstyle_map`` is
    drawn with ``default_markerstyle``), so plotting never fails with KeyError.

    :param config: PlotConfig with the style maps
    :param cells: cell type, e.g. ``"CS"``
    :param media: media, e.g. ``"M9"``
    :param green_intensity: green light intensity, e.g. ``1.4``
    :return: dict with keys ``marker``, ``markercolor``, ``linestyle``, ``linecolor``
    """
    return {
        "marker": config.markerstyle_map.get(cells, config.default_markerstyle),
        "markercolor": config.markercolor_map.get(cells, config.default_markercolor),
        "linestyle": config.linestyle_map.get(media, config.default_linestyle),
        "linecolor": config.line_color_map.get(green_intensity, config.default_linecolor),
    }


def apply_exclusions(df, config, columns=("cells", "media", "green_intensity", "red_intensity")):
    """
    Return a copy of ``df`` without the rows listed in ``config.plot_exclude``.

    Uses the same rule as the old ``if ... in config.plot_exclude[...]``
    checks (exact match). Plot functions filter their data with this once and
    pass the same filtered data to :func:`add_legends`, so the plot and the
    legends always agree.

    :param df: long-format or summary dataframe
    :param config: PlotConfig with a ``plot_exclude`` dict
    :param columns: which exclusion keys to apply. The by-intensity plots only
        apply ``("cells", "media")`` to the means because intensity is the x-axis.
    :return: filtered copy of ``df``
    """
    mask = pd.Series(True, index=df.index)
    for column in columns:
        excluded = config.plot_exclude.get(column, [])
        if column in df.columns and excluded:
            mask &= ~df[column].isin(excluded)
    return df.loc[mask].copy()


def _present_values(df, column):
    """Unique non-NaN values of ``column``; numbers sorted high -> low, text A -> Z."""
    values = pd.Series(df[column].dropna().unique())
    if values.empty:
        return []
    if pd.api.types.is_numeric_dtype(values):
        return sorted(values.tolist(), key=float, reverse=True)
    return sorted(values.tolist(), key=str)


def _handles_for(column, values, config, label_map=None):
    """
    Build legend handles (dummy ``Line2D`` objects) for one column.

    Styles come from :func:`style_for`, so legends match the plot exactly,
    including fallback styles.
    """
    label_map = label_map or {}
    handles = []
    for value in values:
        label = label_map.get(value, f"{value:g}" if isinstance(value, (int, float, np.number)) else str(value))

        if column == "cells":
            style = style_for(config, cells=value)
            handles.append(Line2D([], [], linestyle="None",
                                  marker=style["marker"],
                                  markerfacecolor=style["markercolor"],
                                  markeredgecolor=style["markercolor"],
                                  markersize=6, label=label))
        elif column == "media":
            style = style_for(config, media=value)
            handles.append(Line2D([], [], color="black", linewidth=2,
                                  linestyle=style["linestyle"], label=label))
        elif column == "green_intensity":
            style = style_for(config, green_intensity=value)
            handles.append(Line2D([], [], linewidth=2,
                                  color=style["linecolor"], label=label))
        else:  # any other column: plain grey line
            handles.append(Line2D([], [], color="grey", linewidth=2, label=label))
    return handles


def add_legends(ax, data, config, legends=None, label_maps=None,
                pad=0.02, gap=0.03, skip_single=False, fontsize="small"):
    """
    Add one legend per metadata column, showing only values present in ``data``.

    Legends are placed to the right of the axes, stacked top to bottom. Their
    real sizes are measured after drawing, so they never overlap; if they are
    taller than the axes, they wrap into a second column. The axes are then
    narrowed just enough to make room.

    Call AFTER ``fig.tight_layout()`` (tight_layout would otherwise undo the
    resizing) and save with ``bbox_inches="tight"`` (done by :func:`save_or_show`).

    :param ax: matplotlib Axes that was plotted on
    :param data: the dataframe that was actually plotted (already filtered with
        :func:`apply_exclusions`)
    :param config: PlotConfig with the style maps
    :param legends: dict ``{column: title}`` in display order. Default:
        cells, media, green intensity (:data:`LEGEND_SPECS`).
    :param label_maps: optional ``{column: {value: readable_label}}``,
        e.g. ``{"cells": cell_label_map}``
    :param pad: gap between axes and legends, as a fraction of axes width
    :param gap: vertical gap between legends, as a fraction of axes height
    :param skip_single: if True, omit a legend whose column has only one value
        (e.g. a single media type)
    :param fontsize: legend font size
    :return: list of the Legend objects added
    """
    legends = legends or LEGEND_SPECS
    label_maps = label_maps or {}
    fig = ax.figure

    # --- 1. Create legends with only the present values -----------------------
    created = []
    for column, title in legends.items():
        if column not in data.columns:
            continue
        values = _present_values(data, column)
        if not values or (skip_single and len(values) == 1):
            continue
        handles = _handles_for(column, values, config, label_maps.get(column))
        leg = ax.legend(handles=handles, title=title, loc="upper left",
                        bbox_to_anchor=(1 + pad, 1.0), borderaxespad=0.0,
                        fontsize=fontsize, title_fontsize=fontsize)
        ax.add_artist(leg)   # keep it when the next ax.legend() is called
        # add_artist clips to the axes, which makes bbox_inches="tight" ignore
        # (and crop) legends outside the axes - switch that off
        leg.set_clip_on(False)
        created.append(leg)

    if not created:
        return []

    # --- 2. Measure legend sizes in pixels -------------------------------------
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    sizes = [leg.get_window_extent(renderer) for leg in created]

    # --- 3. Pack legends into columns that fit the axes height -----------------
    axes_h_px = ax.bbox.height
    gap_px = gap * axes_h_px
    columns, current, used = [], [], 0.0
    for leg, box in zip(created, sizes):
        needed = box.height + (gap_px if current else 0)
        if current and used + needed > axes_h_px:   # start a new column
            columns.append(current)
            current, used = [], 0.0
            needed = box.height
        current.append((leg, box))
        used += needed
    columns.append(current)

    col_widths_px = [max(box.width for _, box in col) for col in columns]

    # --- 4. Narrow the axes so all legend columns fit inside the figure --------
    fig_w_px = fig.bbox.width
    pad_px = pad * ax.bbox.width
    total_w_frac = (sum(col_widths_px) + pad_px * (len(columns) + 1)) / fig_w_px
    pos = ax.get_position()
    new_right = min(pos.x1, 1.0 - total_w_frac)
    if new_right > pos.x0 + 0.1:            # never squash the plot below 10 %
        ax.set_position([pos.x0, pos.y0, new_right - pos.x0, pos.height])

    # --- 5. Place each legend (anchors in axes fractions) ----------------------
    ax_w_px, ax_h_px = ax.bbox.width, ax.bbox.height
    x_px = ax_w_px + pad_px
    for col, col_w in zip(columns, col_widths_px):
        y_px = ax_h_px
        for leg, box in col:
            leg.set_bbox_to_anchor((x_px / ax_w_px, y_px / ax_h_px), transform=ax.transAxes)
            y_px -= box.height + gap_px
        x_px += col_w + pad_px

    return created


def save_or_show(fig, title, config, save_image):
    """
    Save the figure to ``<config.save_file_path>/figures/<title>.png`` or show it.

    ``bbox_inches="tight"`` makes sure legends outside the axes are included.
    ``"/"`` in the title is replaced with ``"_div_"`` (not allowed in filenames).
    The figure is always closed afterwards to free memory.

    :param fig: matplotlib Figure
    :param title: plot title, used as the filename
    :param config: PlotConfig (``save_file_path`` is used)
    :param save_image: True -> save PNG; False -> display interactively
    """
    if save_image:
        try:
            folder = Path(config.save_file_path) / "figures"
            folder.mkdir(parents=True, exist_ok=True)
            save_title = title.replace("/", "_div_")
            fig.savefig(folder / f"{save_title}.png", bbox_inches="tight")
        except OSError as error:
            print(f"Could not save image: {error}")
    else:
        plt.show()
    plt.close(fig)


def load_data(filepath, plate_map):
    """
    Load the combined CSV, attach plate metadata, add normalised
    measurements, and summarise replicates.

    Steps:
        1. Read the combined CSV (two-row header: measurement, timestamp) and
           convert timestamps to datetime (UK day-first format).
        2. "Melt" into long format: one row per well x measurement x timepoint,
           tagged with that well's cells / media / green / red intensity.
        3. For every non-OD measurement X, add a derived measurement
           ``"X/OD600"`` = X / OD600 for the same well and timepoint.
           (No blank/background subtraction is applied.)
        4. Average wells with identical metadata to give mean, SD and n.

    Side effect: sets ``DefaultConfig.line_color_map`` to an automatic
    red -> green colour scale based on the green intensities present.

    :param filepath: path to ``*_extracted_combined.csv``
    :param plate_map: DataFrame from :func:`load_plate_metadata`, indexed by well
    :return: ``(sorted_data_df, summary_df)`` - see module docstring for columns.
        ``time`` is hours since the first read of that measurement;
        ``timepoint`` is the 0-based index of the read.
    """

    # Column definitions (documentation only - this empty frame is replaced
    # further down by pd.DataFrame(rows)).
    # TIDY: this block can be removed; "well" is also missing from the list.
    columns = [
                "cells", #str e.g. YX001
                "media", #str e.g. WM-met+
                "green_intensity", #float e.g. 2.8
                "red_intensity", #float e.g. 2.8
                "timestamp", #datetime
                "timepoint", #int, e.g. 1,2. Count of number of times data is collected
                "time", #float in hours
                "measurement", #str e.g. OD600, GFP 395nm
                "value", #float
               ]
    sorted_data_df = pd.DataFrame(columns=columns).astype(object)

    # --- 1. Load data ---------------------------------------------------------
    # header=[0, 1]: first two rows form a (measurement, timestamp) MultiIndex.
    # The third row ("well,,,,") becomes the index name.
    data = pd.read_csv(filepath, header = [0,1], index_col= 0)

    # Rebuild the column MultiIndex so timestamps are real datetimes
    # (needed for time arithmetic). dayfirst=True: "19/11/2026" = 19 Nov.
    timestamps = pd.to_datetime(data.columns.get_level_values("timestamp"), dayfirst=True)
    measurements = data.columns.get_level_values("measurement")
    data.columns = pd.MultiIndex.from_arrays([measurements, timestamps], names=["measurement", "timestamp"])
    rows = []

    # --- 2. Convert to long format -------------------------------------------
    for well in data.index:

        # Metadata for this well from the plate map
        cells = plate_map.loc[well,"cells"]
        media = plate_map.loc[well,"media"]
        green_intensity = plate_map.loc[well,"green_intensity"]
        red_intensity = plate_map.loc[well,"red_intensity"]

        # One row per (measurement, timestamp) column
        for measurement, timestamp in data.columns:
            timestamps = data[measurement].columns       # all read times for this measurement
            timepoint = timestamps.get_loc(timestamp)    # 0-based read number
            initial_time = timestamps.min()              # first read = time zero
            time = (timestamp - initial_time).total_seconds() /3600   # hours
            value = data.loc[well, (measurement, timestamp)]
            rows.append({
                "well": well,
                "cells": cells,
                "media": media,
                "green_intensity": green_intensity,
                "red_intensity": red_intensity,
                "timestamp": timestamp,
                "timepoint": timepoint,
                "time" : time,
                "measurement": measurement,
                "value": value,
            })

    sorted_data_df = pd.DataFrame(rows)

    # --- 3. Normalise each signal by OD600 ------------------------------------
    # e.g. "GFP 395nm" -> new measurement "GFP 395nm/OD600"
    # TIDY: row-by-row lookup is slow (O(n^2)); a merge on (well, timepoint)
    # would do the same thing in one step - see suggestions.
    select_od_data = sorted_data_df.loc[sorted_data_df["measurement"] == "OD600"].copy()
    new_rows = []
    for measurement in sorted_data_df["measurement"].unique():
        if measurement == "OD600":
            continue
        else:
            select_measurement_data = sorted_data_df.loc[sorted_data_df["measurement"] == measurement].copy()
            calculated_measurement = f"{measurement}/OD600"
            for index, measure_row in select_measurement_data.iterrows():
                well = measure_row["well"]
                timepoint = measure_row["timepoint"]

                # OD600 reading for the same well and read number
                od_row = select_od_data[(select_od_data["well"] == well) & (select_od_data["timepoint"] == timepoint)]
                od_value = od_row["value"].iloc[0]

                # Fluorescence / OD (raw values, no blank subtraction)
                measure_value = measure_row["value"]
                ratio = measure_value / od_value

                # Copy the row (keeps metadata/time) and replace measurement + value
                new_row = measure_row.copy()
                new_row["measurement"] = calculated_measurement
                new_row["value"] = ratio
                new_rows.append(new_row)

    calculated_data_df = pd.DataFrame(new_rows)
    sorted_data_df = pd.concat([sorted_data_df, calculated_data_df], ignore_index= True)

    # --- 4. Average replicate wells -------------------------------------------
    # Wells with the same cells / media / green / red at the same read are
    # treated as replicates. std is the sample SD (ddof=1): NaN if n = 1.
    # Note: groupby drops rows where any key is NaN (e.g. empty wells with no
    # metadata).
    # TIDY: grouping by "timestamp" makes "min_timestamp" identical to it.
    summary_df = sorted_data_df.groupby(["cells","media","green_intensity","red_intensity","timestamp","time","timepoint","measurement"], as_index=False).agg(
        min_timestamp = ("timestamp", "min"),
        mean = ("value","mean"),
        std = ("value","std"),
        count = ("value","count")
    )

    # Automatic colour per green intensity (overridden in main() if the user
    # supplies their own line_color_map)
    DefaultConfig.line_color_map = build_intensity_colour_map(sorted_data_df["green_intensity"].unique())

    return sorted_data_df, summary_df


def plot_timecourse(raw_dataframe, summary_dataframe, measurement, plot_type, ylabel: str | None = None, title: str | None = None, title_extra: str = "",
                            xlabel = "Time (hrs)",
                            config: PlotConfig = DefaultConfig, save_image = False,
                           ):
    """
    Plot one measurement against time, one line per well or condition.

    Line colour = green intensity, marker = cell type, line style = media.
    Legends (cell type, media, green intensity) are drawn to the right and
    list only the conditions actually plotted.

    :param raw_dataframe: long-format data (``sorted_data_df`` from :func:`load_data`)
    :param summary_dataframe: replicate means/SDs (``summary_df`` from :func:`load_data`)
    :param measurement: measurement to plot, e.g. ``"OD600"`` or ``"GFP 395nm/OD600"``
    :param plot_type: ``"average"`` - mean +/- SD per condition;
        ``"all"`` - every well individually;
        ``"both"`` - wells faded (alpha 0.2) under the means
    :param ylabel: y-axis label (default: ``measurement``)
    :type ylabel: str | None
    :param title: title override; ``title_extra`` is still appended
    :type title: str | None
    :param title_extra: text appended to the end of the title (and filename)
    :type title_extra: str
    :param xlabel: x-axis label
    :param config: styling/output settings
    :type config: PlotConfig
    :param save_image: True -> save PNG to ``<save_file_path>/figures/<title>.png``;
        False -> display interactively
    """

    if ylabel is None:
        ylabel = measurement

    if title is None:
        title = f"{ylabel} - {plot_type} - {title_extra}"
    else:
        title = title + " - " + title_extra

    # Keep only the requested measurement, then remove excluded conditions.
    # The same filtered data is used for plotting AND for the legends.
    summary_data = summary_dataframe.loc[summary_dataframe["measurement"] == measurement].copy()
    raw_data = raw_dataframe.loc[raw_dataframe["measurement"] == measurement].copy()
    summary_data = apply_exclusions(summary_data, config)
    raw_data = apply_exclusions(raw_data, config)

    plotted = []  # dataframes actually drawn -> used to build the legends

    fig, axs = plt.subplots()

    # --- Individual wells -----------------------------------------------------
    if plot_type == "all" or plot_type == "both":
        # Fade individual wells when means are drawn on top
        if plot_type == "both":
            select_alpha = 0.2
        else:
            select_alpha = 1.0

        for well, group in raw_data.groupby("well"):
            # Metadata is the same on every row of a well, so read from the first
            row = group.iloc[0]
            style = style_for(config, row["cells"], row["media"], row["green_intensity"])

            axs.plot(group["time"], group["value"],
                        color = style["linecolor"],
                        marker = style["marker"],
                        markerfacecolor = style["markercolor"],
                        markeredgecolor = style["markercolor"], markersize = 3.0,
                        linestyle = style["linestyle"], linewidth = 1.0,
                        alpha = select_alpha,
                        )
        plotted.append(raw_data)

    # --- Condition means +/- SD ------------------------------------------------
    if plot_type == "average" or plot_type == "both":
        for (cells, media, green_intensity, red_intensity), group in summary_data.groupby(["cells","media","green_intensity","red_intensity"]):
            style = style_for(config, cells, media, green_intensity)

            axs.errorbar(group["time"], group["mean"],
                        yerr = group["std"], capsize = 2.0,
                        color = style["linecolor"],
                        marker = style["marker"],
                        markerfacecolor = style["markercolor"],
                        markeredgecolor = style["markercolor"], markersize = 3.0,
                        linestyle = style["linestyle"], linewidth = 1.0,
                        alpha = 1.0,
                        )
        plotted.append(summary_data)

    # --- Labels, legends, save -------------------------------------------------
    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_title(title)
    fig.tight_layout()   # must come BEFORE add_legends
    if plotted:
        add_legends(axs, pd.concat(plotted, ignore_index=True), config)

    save_or_show(fig, title, config, save_image)


def plot_by_intensity(raw_dataframe, summary_dataframe, measurement, plot_type, timepoints: int | list[int], ylabel: str | None = None, title: str | None = None,
                              title_extra: str = "",
                              xlabel = "Green light intensity %",
                                config: PlotConfig = DefaultConfig, save_image = False):
    """
    Dose-response plot: one measurement against green light intensity at a
    chosen read (or reads).

    Intensity is shown as a percentage of 2.8 (assumed maximum intensity).
    The condition green = 2.8 & red = 0 is shifted 10 % to the right of the
    other points so it is visually separated.
    Points and connecting line segments are coloured red -> green by intensity.
    Legends (cell type, media) list only the conditions actually plotted.

    :param raw_dataframe: long-format data (``sorted_data_df``)
    :param summary_dataframe: replicate means/SDs (``summary_df``)
    :param measurement: measurement to plot, e.g. ``"GFP 395nm/OD600"``
    :param plot_type: ``"average"``, ``"all"`` or ``"both"`` (see :func:`plot_timecourse`)
    :param timepoints: 0-based read index (0 = first read), or a list of them.
        A single int uses the red -> green colormap; a list switches to viridis
        (multi-timepoint colouring is not yet implemented, see TODO).
    :param ylabel: y-axis label (default: ``measurement``)
    :type ylabel: str | None
    :param title: title override; ``title_extra`` is still appended
    :type title: str | None
    :param title_extra: text appended to the end of the title (and filename)
    :type title_extra: str
    :param xlabel: x-axis label (drawn in green)
    :param config: styling/output settings
    :type config: PlotConfig
    :param save_image: True -> save PNG; False -> display interactively
    """
    if ylabel is None:
        ylabel = measurement

    if title is None:
        title = f"{ylabel} - by intensity - {plot_type} - {title_extra}"
    else:
        title = title + " - " + title_extra

    # Choose colormap depending on single vs multiple timepoints
    if isinstance(timepoints, int):  # single timepoint allows gradient colours to be used
        timepoints = [timepoints]
        colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
        cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)
        timepoints_is_array = False
    else:
        cmap = plt.get_cmap("viridis")  # Use a different colormap for multiple timepoints
        timepoints_is_array = True
        #TODO - implement different line colours for multiple timepoints

    plotted = []  # dataframes actually drawn -> used to build the legends

    fig, axs = plt.subplots()
    for timepoint in timepoints:
        # Keep only this measurement at this read
        raw_data = raw_dataframe.loc[(raw_dataframe["measurement"] == measurement) & (raw_dataframe["timepoint"] == timepoint)].copy()
        summary_data = summary_dataframe.loc[(summary_dataframe["measurement"] == measurement) & (summary_dataframe["timepoint"] == timepoint)].copy()

        # Remove excluded conditions. Individual wells: all four keys (as before).
        # Means: only cells/media, because intensity is the x-axis.
        raw_data = apply_exclusions(raw_data, config)
        summary_data = apply_exclusions(summary_data, config, columns=("cells", "media"))

        # Convert intensities to % of maximum (2.8 = 100 %)
        # TIDY: 2.8 is hard-coded in several places - make it a constant.
        raw_data["green_intensity_percentage"] = raw_data["green_intensity"].apply(lambda x: 100*x/2.8)
        raw_data["red_intensity_percentage"] = raw_data["red_intensity"].apply(lambda x: 100*x/2.8)
        summary_data["green_intensity_percentage"] = summary_data["green_intensity"].apply(lambda x: 100*x/2.8)
        summary_data["red_intensity_percentage"] = summary_data["red_intensity"].apply(lambda x: 100*x/2.8)

        select_alpha = 1.0

        # --- Individual wells (scatter only, not joined) ---------------------
        if plot_type == "all" or plot_type == "both":
            if plot_type == "both":
                select_alpha = 0.2

            for well, group in raw_data.groupby("well"):
                row = group.iloc[0]
                style = style_for(config, cells=row["cells"])

                x = group["green_intensity_percentage"].values.copy()
                y = group["value"].values.copy()

                # Move the "full green, no red" point (green=2.8 & red=0) to the right
                is_final = ((group["green_intensity"] == 2.8) &(group["red_intensity"] == 0))
                x[is_final.values] = x.max() + 10

                # Sort by x so points are in intensity order
                order = np.argsort(x)
                x = x[order]
                y = y[order]

                # No line segments: individual wells are not connected
                norm = plt.Normalize(0, 100)
                # TIDY: passing both `c=` (colormap) and `facecolor=` is
                # contradictory - the fixed facecolor may override the gradient.
                scatter = axs.scatter(x, y, c=row["green_intensity_percentage"], cmap=cmap, norm=norm, s=20, zorder=5,
                                    marker = style["marker"],
                                    facecolor = style["markercolor"],
                                    edgecolor = style["markercolor"],
                                    alpha = select_alpha,
                                    )
            plotted.append(raw_data)

        # --- Condition means: gradient line + points + error bars -------------
        if plot_type == "average" or plot_type == "both":
                for (cells, media), group in summary_data.groupby(["cells","media"]):
                    style = style_for(config, cells=cells, media=media)

                    x = group["green_intensity_percentage"].values.copy()
                    y = group["mean"].values.copy()
                    yerr = group["std"].values.copy()

                    # Move the "full green, no red" point to the right
                    is_final = ((group["green_intensity"] == 2.8) &(group["red_intensity"] == 0))
                    x[is_final.values] = x.max() + 10

                    order = np.argsort(x)
                    x = x[order]
                    y = y[order]
                    yerr = yerr[order]

                    # Build consecutive (x, y) -> (x, y) segments so each segment
                    # can take its own colour from the colormap
                    points = np.array([x, y]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)

                    # Map x range onto the colormap (low = red, high = green)
                    norm = plt.Normalize(x.min(), x.max())
                    lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
                    lc.set_linestyle(style["linestyle"])
                    lc.set_alpha(select_alpha)
                    lc.set_array(x)
                    axs.add_collection(lc)

                    # Points on top of the line
                    scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5,
                                        marker = style["marker"],
                                        facecolor = style["markercolor"],
                                        edgecolor = style["markercolor"],
                                        alpha = select_alpha,
                                        )

                    # Error bars only (fmt='none' = no extra markers/lines)
                    axs.errorbar(x, y, yerr=yerr, fmt='none',
                                ecolor= style["markercolor"], alpha=1.0, capsize=2.0,
                                )
                plotted.append(summary_data)

    # --- Labels, legends (no intensity legend - it's the x-axis), save -------
    x_axis = axs.set_xlabel(xlabel)
    x_axis.set_color("green")
    axs.set_ylabel(ylabel)
    axs.set_title(title)

    fig.tight_layout()   # must come BEFORE add_legends
    if plotted:
        add_legends(axs, pd.concat(plotted, ignore_index=True), config,
                    legends={"cells": "Cell type", "media": "Media"})

    save_or_show(fig, title, config, save_image)


def plot_by_intensity_foldchange(summary_dataframe, measurement, timepoints: int | list[int], ylabel: str | None = None, title: str | None = None,
                              title_extra: str = "",
                              xlabel = "Green light intensity %",
                                config: PlotConfig = DefaultConfig, save_image = False):
    """
    Dose-response plot of fold change relative to the red-light control.

    For each cell type / media combination::

        fold change = mean(green intensity X) / mean(green = 0, red = 2.8)

    with the SD propagated in quadrature
    (relative errors added: ``sqrt((sd_g/g)^2 + (sd_r/r)^2) * fold``).
    Only conditions with green > 0 are plotted. Styling as in
    :func:`plot_by_intensity`; legends list only the plotted cell types/media.

    Requires a "100 % red" control (green = 0, red = 2.8) for every
    cell/media group; otherwise the division fails.

    :param summary_dataframe: replicate means/SDs (``summary_df``)
    :param measurement: measurement to plot, e.g. ``"GFP 395nm/OD600"``
    :param timepoints: 0-based read index, or a list of them
    :param ylabel: y-axis label (default: ``measurement``)
    :type ylabel: str | None
    :param title: title override; ``title_extra`` is still appended
    :type title: str | None
    :param title_extra: text appended to the end of the title (and filename)
    :type title_extra: str
    :param xlabel: x-axis label (drawn in green)
    :param config: styling/output settings
    :type config: PlotConfig
    :param save_image: True -> save PNG; False -> display interactively
    """
    if ylabel is None:
        ylabel = measurement

    if title is None:
        title = f"{ylabel} - by intensity - average - {title_extra}"
    else:
        title = title + " - " + title_extra

    if isinstance(timepoints, int):  # single timepoint allows gradient colours to be used
        timepoints = [timepoints]
        colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
        cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)
        timepoints_is_array = False
    else:
        cmap = plt.get_cmap("viridis")  # Use a different colormap for multiple timepoints
        timepoints_is_array = True
        #TODO - implement different line colours for multiple timepoints

    plotted = []  # dataframes actually drawn -> used to build the legends

    fig, axs = plt.subplots()
    for timepoint in timepoints:
        # Keep only this measurement at this read, drop excluded cells/media,
        # add % intensity columns
        summary_data = summary_dataframe.loc[(summary_dataframe["measurement"] == measurement) & (summary_dataframe["timepoint"] == timepoint)].copy()
        summary_data = apply_exclusions(summary_data, config, columns=("cells", "media"))
        summary_data["green_intensity_percentage"] = summary_data["green_intensity"].apply(lambda x: 100*x/2.8)
        summary_data["red_intensity_percentage"] = summary_data["red_intensity"].apply(lambda x: 100*x/2.8)

        for (cells, media), group in summary_data.groupby(["cells","media"]):
            style = style_for(config, cells=cells, media=media)

            # Numerator: every condition with some green light
            green_data = group.loc[group["green_intensity"] > 0].copy()
            # .copy(): newer pandas returns a read-only view, and x is edited below
            x = green_data["green_intensity_percentage"].values.copy()

            y_green = green_data["mean"].values
            yerr_green = green_data["std"].values

            # Denominator: the red-only control (green = 0, red = 2.8)
            y_red = group.loc[(group["red_intensity"] == 2.8) & (group["green_intensity"] == 0), "mean"].values.copy()
            yerr_red = group.loc[(group["red_intensity"] == 2.8) & (group["green_intensity"] == 0), "std"].values.copy()

            y = y_green / y_red  # fold change
            yerr = np.sqrt((yerr_green / y_green)**2 + (yerr_red / y_red)**2) * y  # propagated error

            # Move the "full green, no red" point to the right
            is_final = ((green_data["green_intensity"] == 2.8) &(green_data["red_intensity"] == 0))
            x[is_final.values] = x.max() + 10

            order = np.argsort(x)
            x = x[order]
            y = y[order]
            yerr = yerr[order]

            # Gradient-coloured connecting line
            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            norm = plt.Normalize(x.min(), x.max())
            lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0, alpha=1.0)
            lc.set_linestyle(style["linestyle"])
            lc.set_array(x)
            axs.add_collection(lc)

            # Points
            scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5,
                                marker = style["marker"],
                                facecolor = style["markercolor"],
                                edgecolor = style["markercolor"],
                                alpha = 1.0,
                                )

            # Error bars
            axs.errorbar(x, y, yerr=yerr, fmt='none',
                        ecolor= style["markercolor"], alpha=1.0, capsize=2.0,
                        )
        plotted.append(summary_data)

    # --- Labels, legends, save -------------------------------------------------
    x_axis = axs.set_xlabel(xlabel)
    x_axis.set_color("green")
    axs.set_ylabel(ylabel)
    axs.set_title(title)

    fig.tight_layout()   # must come BEFORE add_legends
    if plotted:
        add_legends(axs, pd.concat(plotted, ignore_index=True), config,
                    legends={"cells": "Cell type", "media": "Media"})

    save_or_show(fig, title, config, save_image)


def plot_timecourse_custom(dataframe, y_data, plot_type, ylabel: str | None = None, title: str | None = None, title_extra: str = "",
                            xlabel = "Time (hrs)",
                            config: PlotConfig = DefaultConfig, save_image = False,
                            row_filter = None,
                           ):
    """
    LEGACY timecourse plot for the older *wide* data format.

    Expects one row per condition with list-valued columns
    ``f"{y_data}_timepoints"``, ``f"{y_data}_average"``, ``f"{y_data}_std"``
    and ``f"{y_data}_raw_array"``. This format is no longer produced by
    :func:`load_data`. (The former NameError in the legend code is fixed by
    using :func:`add_legends`.)

    :param dataframe: wide-format data, one row per condition
    :param y_data: base column name, e.g. ``"OD600"``
    :param plot_type: ``"average"``, ``"all"``, or ``"both"``
        (``"both"`` uses hard-coded colours for JBL137 vs other strains, so no
        green-intensity legend is drawn in that mode)
    :param ylabel: y-axis label (default: ``y_data``)
    :type ylabel: str | None
    :param title: title override; ``title_extra`` is still appended
    :type title: str | None
    :param title_extra: text appended to the end of the title
    :type title_extra: str
    :param xlabel: x-axis label
    :param config: styling/output settings
    :type config: PlotConfig
    :param save_image: True -> save PNG; False -> display interactively
    :param row_filter: optional function ``row -> bool``; rows returning False
        are skipped
    """

    if ylabel is None:
        ylabel = y_data

    if title is None:
        title = f"{ylabel} - {plot_type} - {title_extra}"
    else:
        title = title + " - " + title_extra

    # Remove excluded conditions and apply the optional user filter up front,
    # so the legends can be built from exactly the rows that are plotted
    dataframe = apply_exclusions(dataframe, config)
    if row_filter is not None and not dataframe.empty:
        dataframe = dataframe[dataframe.apply(row_filter, axis=1)]

    fig, axs = plt.subplots()
    for row_index, row in dataframe.iterrows():

        index_name_cells = row["cells"]
        index_name_media = row["media"]
        index_name_green_intensity = row["green_intensity"]
        style = style_for(config, index_name_cells, index_name_media, index_name_green_intensity)

        # Optional per-intensity transparency
        if config.alpha_used is True:
            alpha = config.line_alpha_map.get(index_name_green_intensity, config.default_linealpha)
        else:
            alpha = 1

        if plot_type == "average":
            axs.errorbar(row[f"{y_data}_timepoints"], row[f"{y_data}_average"],
                        yerr = row[f"{y_data}_std"], capsize = 2.0,
                        color = style["linecolor"],
                        marker = style["marker"],
                        markerfacecolor = style["markercolor"],
                        markeredgecolor = style["markercolor"], markersize = 3.0,
                        linestyle = style["linestyle"], linewidth = 1.0,
                        alpha = alpha,
                        )

        elif plot_type == "all":
            # One line per replicate
            for repeat in row[f"{y_data}_raw_array"]:
                axs.plot(row[f"{y_data}_timepoints"], repeat,
                            color = style["linecolor"],
                            marker = style["marker"],
                            markerfacecolor = style["markercolor"],
                            markeredgecolor = style["markercolor"], markersize = 3.0,
                            linestyle = style["linestyle"], linewidth = 1.0,
                            alpha = alpha,
                            )

        elif plot_type == "both":

            # Hard-coded colours from an earlier JBL137 experiment
            if row["cells"] == "JBL137":
                custom_line_color = {0.28: "blue", 0.0: "red", 2.8: "blue", 0.028:"blue"}
            else:
                custom_line_color = {0.28: "black", 0.0: "black", 2.8: "black", 0.028:"blue"}

            # Mean +/- SD
            axs.errorbar(row[f"{y_data}_timepoints"], row[f"{y_data}_average"],
                        yerr = row[f"{y_data}_std"], capsize = 2.0,
                        color = custom_line_color[index_name_green_intensity],
                        marker = style["marker"],
                        markerfacecolor = custom_line_color[index_name_green_intensity],
                        markeredgecolor = custom_line_color[index_name_green_intensity], markersize = 3.0,
                        linestyle = style["linestyle"], linewidth = 1.0,
                        alpha = 1.0,
                        )

            # Faded individual replicates (skipped for JBL001 and media)
            for repeat in row[f"{y_data}_raw_array"]:
                if row["cells"] == "JBL001" or row["cells"] == "media":
                    continue
                axs.plot(row[f"{y_data}_timepoints"], repeat,
                            color = custom_line_color[index_name_green_intensity],
                            linestyle = style["linestyle"], linewidth = 1.0,
                            alpha = 0.2,
                            )

    # --- Labels, legends, save -------------------------------------------------
    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_title(title)
    fig.tight_layout()   # must come BEFORE add_legends

    if plot_type == "both":
        # colours are hard-coded above, so an intensity legend would be misleading
        add_legends(axs, dataframe, config, legends={"cells": "Cell type", "media": "Media"})
    else:
        add_legends(axs, dataframe, config)

    save_or_show(fig, title, config, save_image)


def load_plate_metadata(filepath):
    """
    Read the experiment plate map (which conditions are in which well).

    The Excel sheet uses the same layout as the downstairs plate reader export:
    a ``"Label: <field>"`` row, then 8 rows starting with A-H, each followed by
    12 values (columns 1-12). Repeat for every field, e.g.::

        Label: cells
        A | CS | CS | KS | ...
        ...
        H | ...
        Label: green_intensity
        A | 0  | 0.028 | ...

    Expected fields: ``cells``, ``media``, ``green_intensity``,
    ``red_intensity``. Any other labels are added as extra columns.

    :param filepath: path to ``experiment_input.xlsx``
    :return: DataFrame indexed by well (``A1``-``H12``) with columns
        ``cells``, ``media``, ``green_intensity``, ``red_intensity``,
        ``long_name`` (= ``"<cells>_<media>_<green>_<red>"``)
    """
    # TIDY: these plate definitions duplicate those in plate_reader_extraction.py
    key_rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
    key_columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    key_wells = [str(row)+str(col) for row in key_rows for col in key_columns]

    # Empty 96-well table to fill in
    plate_table_df = pd.DataFrame(index = key_wells, columns = ["cells", "media", "green_intensity", "red_intensity", "long_name"])
    plate_table_df.index.name = "well"

    # Scan the sheet row by row
    data = pd.read_excel(filepath, 0, header=None, engine="calamine")
    for row_num, row in enumerate(data.values):
        first_cell = str(row[0]).strip()

        # "Label: cells" -> current field = "cells"
        if "Label:" in first_cell:
            data_name = row[0].split(":")[-1].strip()

        # Plate row A-H: write 12 values into the current field's column
        if first_cell in key_rows:
            for col_index, value in enumerate(row[1:13]):
                well_id = first_cell + str(col_index + 1)
                plate_table_df.at[well_id, data_name] = value

    # Combined descriptive name, e.g. "CS_M9_1.4_0"
    for well in plate_table_df.index:
        cell = plate_table_df.at[well, "cells"]
        media = plate_table_df.at[well, "media"]
        green_intensity = plate_table_df.at[well, "green_intensity"]
        red_intensity = plate_table_df.at[well, "red_intensity"]
        long_name = f"{cell}_{media}_{green_intensity}_{red_intensity}"
        plate_table_df.at[well, "long_name"] = long_name

    return plate_table_df


def print_style_templates(plate_df):
    """
    Print the unique metadata values in the plate map, plus copy-paste-ready
    Python dictionaries for the USER VARIABLES block.

    Printed templates: ``cell_label_map``, ``markerstyle_map``,
    ``markercolor_map``, ``linestyle_map``, ``line_color_map`` (evenly spaced
    red -> green gradient by rank of green intensity) and ``plot_exclude``
    (all values commented out).

    :param plate_df: DataFrame from :func:`load_plate_metadata`
    """

    # Unique values of each metadata field
    cells = sorted(plate_df["cells"].dropna().unique())
    media = sorted(plate_df["media"].dropna().unique())
    green = sorted(
        pd.to_numeric(
            plate_df["green_intensity"].dropna()
        ).unique()
    )
    red = sorted(
        pd.to_numeric(
            plate_df["red_intensity"].dropna()
        ).unique()
    )

    print("\n### UNIQUE PLATE METADATA ###")
    print("Cells:", cells)
    print("Media:", media)
    print("Green intensities:", green)
    print("Red intensities:", red)

    print("\n### COPY AND PASTE STYLE MAPS ###\n")

    # Human-readable cell labels (defaults to the code itself)
    print("\n### HUMAN-READABLE CELL LABELS ###\n")
    print("cell_label_map = {")

    cells = sorted(plate_df["cells"].dropna().unique())  # TIDY: duplicate of above

    for cell in cells:
        print(f"    {repr(cell)}: {repr(str(cell))},")

    print("}\n")

    def print_dict(name, keys, default):
        """Print ``name = {key: default, ...}`` as valid Python source."""
        print(f"{name} = {{")
        for key in keys:
            print(f"    {repr(key)}: {repr(default)},")
        print("}\n")

    print_dict("markerstyle_map", cells, "o")
    print_dict("markercolor_map", cells, "black")
    print_dict("linestyle_map", media, "solid")

    # Line colours: evenly spaced by RANK (not value) from red (lowest) to green (highest)
    print("\nline_color_map = {")

    green_values = sorted(
        pd.to_numeric(
            plate_df["green_intensity"].dropna()
        ).unique()
    )  # TIDY: same as `green` above

    min_value = min(green_values)
    max_value = max(green_values)

    for value in green_values:

        if max_value == min_value:
            fraction = 0.5
        else:
            fraction = list(green_values).index(value) / (len(green_values) - 1)

        # Linear gradient: red -> green
        # TIDY: `red` / `green` here overwrite the lists of the same name above
        red = 1.0 - fraction
        green = fraction
        blue = 0.0

        print(
            f"    {float(value)!r}: "
            f"({red:.3f}, {green:.3f}, {blue:.3f}),"
        )

    print("}\n")

    # Exclusion template: every value listed but commented out
    print("\nplot_exclude = {")

    for column in [
        "cells",
        "media",
        "green_intensity",
        "red_intensity"
    ]:
        values = sorted(plate_df[column].dropna().unique())

        print(f'    "{column}": [')

        for value in values:
            # Convert NumPy numbers to standard Python numbers for clean printing
            if isinstance(value, np.generic):
                value = value.item()

            print(f"        # {repr(value)},")

        print("    ],")

    print("}")


### BAR CHARTS AND STATISTICS ###
#
# Workflow:
#   _replicate_values()        one value per biological replicate
#   compare_conditions_stats() CS vs KS test at every (timepoint, intensity)
#   plot_comparison_bars()     grouped bars (raw or fold change) + p-values
#
# Replicates: wells with a colony suffix (CS_1, CS_2 ...) are averaged per
# colony first and each colony counts as one replicate. Wells without a suffix
# (e.g. plain "CS") each count as one replicate. See `replicate` below.

def timepoint_at_hours(raw_dataframe, hours):
    """
    Return the 0-based read index (timepoint) closest to a time in hours.

    Example: ``timepoint_at_hours(sorted_data_df, 12)`` -> 4 for this plate
    (reads at ~0, 3.6, 6, 8.2, 12.1, 15, 24.1 h).

    :param raw_dataframe: long-format data (``sorted_data_df``)
    :param hours: time in hours
    :return: int timepoint index
    """
    times = raw_dataframe.groupby("timepoint")["time"].mean()
    return int((times - hours).abs().idxmin())


def _replicate_values(raw_dataframe, measurement, config=DefaultConfig, replicate="auto"):
    """
    One value per biological replicate for every timepoint x intensity x cell type.

    * Exclusions in ``config.plot_exclude`` for media / green / red intensity
      are applied. For cells, only individual colonies (names ending in
      ``_<number>``, e.g. ``"CS_3"``) are excluded here - which cell types are
      shown is chosen with ``conditions`` in :func:`plot_comparison_bars`, and
      WT / media must stay available as fold-change references.
    * ``cell_type`` = cell name without colony suffix (``CS_1`` -> ``CS``).

    :param raw_dataframe: long-format data (``sorted_data_df``)
    :param measurement: e.g. ``"GFP 395nm/OD600"`` or ``"OD600"``
    :param config: PlotConfig (uses ``plot_exclude``)
    :param replicate: what counts as one replicate:
        ``"auto"``   - colony if the cell name has a ``_<n>`` suffix, otherwise well;
        ``"well"``   - every well;
        ``"colony"`` - every distinct cell name (wells averaged per name)
    :return: DataFrame with columns ``timepoint, time, green_intensity,
        cell_type, replicate_id, value``
    """
    df = raw_dataframe.loc[raw_dataframe["measurement"] == measurement].copy()
    if df.empty:
        raise ValueError(f"No data for measurement '{measurement}'")

    df = apply_exclusions(df, config, columns=("media", "green_intensity", "red_intensity"))

    df["green_intensity"] = pd.to_numeric(df["green_intensity"])
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df["colony"] = df["cells"].astype(str)
    df["cell_type"] = df["colony"].str.replace(r"_\d+$", "", regex=True)

    # Exclude individual colonies only (e.g. "CS_3")
    excluded_colonies = [c for c in config.plot_exclude.get("cells", [])
                         if re.search(r"_\d+$", str(c))]
    df = df[~df["colony"].isin(excluded_colonies)]

    if df["media"].nunique() > 1:
        warnings.warn("More than one media type present - they are pooled. "
                      "Exclude all but one with plot_exclude['media'].")

    # Choose the replicate unit
    if replicate == "well":
        df["replicate_id"] = df["well"]
    elif replicate == "colony":
        df["replicate_id"] = df["colony"]
    elif replicate == "auto":
        has_suffix = df["colony"] != df["cell_type"]
        df["replicate_id"] = np.where(has_suffix, df["colony"], df["well"])
    else:
        raise ValueError("replicate must be 'auto', 'well' or 'colony'")

    # Average wells belonging to the same replicate (technical replicates)
    reps = (df.groupby(["timepoint", "green_intensity", "cell_type", "replicate_id"], as_index=False)
              .agg(time=("time", "mean"), value=("value", "mean")))
    return reps


def _run_test(a, b, test):
    """
    Two-sample test between arrays ``a`` and ``b``. Returns ``(statistic, p)``,
    or ``(nan, nan)`` when there are too few values (n < 2 in either group).

    :param test: ``"welch"`` (t-test, unequal variances - default),
        ``"student"`` (t-test, equal variances),
        ``"welch_log"`` (Welch t-test on log10 values - for ratio/fold data;
        requires all values > 0),
        ``"mannwhitney"`` (rank-based, no normality assumption)
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan

    if test == "welch":
        result = sps.ttest_ind(a, b, equal_var=False)
    elif test == "student":
        result = sps.ttest_ind(a, b, equal_var=True)
    elif test == "welch_log":
        if (a <= 0).any() or (b <= 0).any():
            warnings.warn("welch_log needs positive values - returning NaN")
            return np.nan, np.nan
        result = sps.ttest_ind(np.log10(a), np.log10(b), equal_var=False)
    elif test == "mannwhitney":
        result = sps.mannwhitneyu(a, b, alternative="two-sided")
    else:
        raise ValueError("test must be 'welch', 'student', 'welch_log' or 'mannwhitney'")
    return float(result.statistic), float(result.pvalue)


def _adjust_pvalues(p, method="holm"):
    """
    Multiple-comparison correction. NaNs are ignored and kept as NaN.

    :param p: array of p-values
    :param method: ``"holm"`` (controls family-wise error, recommended for a
        handful of planned comparisons), ``"bh"`` (Benjamini-Hochberg false
        discovery rate, for many comparisons) or ``None`` (no correction)
    :return: array of adjusted p-values
    """
    p = np.asarray(p, dtype=float)
    adjusted = np.full_like(p, np.nan)
    valid = ~np.isnan(p)
    pv = p[valid]
    m = len(pv)
    if m == 0 or method is None:
        return p.copy()

    order = np.argsort(pv)
    ranked = pv[order]
    if method == "holm":
        adj = np.maximum.accumulate((m - np.arange(m)) * ranked)
    elif method == "bh":
        adj = np.minimum.accumulate((m / np.arange(m, 0, -1)) * ranked[::-1])[::-1]
    else:
        raise ValueError("method must be 'holm', 'bh' or None")
    adj = np.minimum(adj, 1.0)

    out = np.empty(m)
    out[order] = adj
    adjusted[valid] = out
    return adjusted


def _p_to_stars(p):
    """``*** p<0.001, ** p<0.01, * p<0.05, ns otherwise`` (blank if p is NaN)."""
    if np.isnan(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def compare_conditions_stats(raw_dataframe, measurement, group_a="CS", group_b="KS",
                             timepoints=None, intensities=None,
                             test="welch", correction="holm", replicate="auto",
                             config=DefaultConfig):
    """
    Test ``group_a`` vs ``group_b`` at every timepoint x green intensity.

    One test per (timepoint, intensity) on the replicate values from
    :func:`_replicate_values`. P-values are then corrected for multiple
    comparisons across ALL rows returned (so restrict ``timepoints`` /
    ``intensities`` to what you actually want to compare).

    Note: dividing both groups by the same WT or media mean does not change
    the test, so these p-values also apply to fold-change bars.

    :param raw_dataframe: long-format data (``sorted_data_df``)
    :param measurement: e.g. ``"GFP 395nm/OD600"``
    :param group_a: first cell type (default ``"CS"``)
    :param group_b: second cell type (default ``"KS"``)
    :param timepoints: int / list of 0-based read indices; None = all
    :param intensities: float / list of green intensities; None = all
    :param test: see :func:`_run_test` (default ``"welch"``), or ``"anova"``
        to use :func:`anova_stats` (two-way ANOVA across intensities; the
        omnibus table is attached as ``result.attrs["anova"]``)
    :param correction: ``"holm"`` (default), ``"bh"`` or ``None``
    :param replicate: ``"auto"``, ``"well"`` or ``"colony"`` (see :func:`_replicate_values`)
    :param config: PlotConfig (exclusions)
    :return: DataFrame, one row per timepoint x intensity, columns
        ``timepoint, time, green_intensity, n_a, n_b, mean_a, mean_b,
        sd_a, sd_b, ratio_a_over_b, statistic, p, p_adj, stars, test``
    """
    if test == "anova":
        posthoc, omnibus = anova_stats(raw_dataframe, measurement, group_a, group_b,
                                       across="intensity", timepoints=timepoints,
                                       intensities=intensities, correction=correction,
                                       replicate=replicate, config=config)
        posthoc.attrs["anova"] = omnibus
        return posthoc

    reps = _replicate_values(raw_dataframe, measurement, config, replicate)

    if timepoints is not None:
        timepoints = [timepoints] if np.isscalar(timepoints) else list(timepoints)
        reps = reps[reps["timepoint"].isin(timepoints)]
    if intensities is not None:
        intensities = [intensities] if np.isscalar(intensities) else list(intensities)
        reps = reps[np.isclose(reps["green_intensity"].to_numpy()[:, None], intensities).any(axis=1)]

    rows = []
    for (timepoint, intensity), group in reps.groupby(["timepoint", "green_intensity"]):
        a = group.loc[group["cell_type"] == group_a, "value"].to_numpy()
        b = group.loc[group["cell_type"] == group_b, "value"].to_numpy()
        if len(a) == 0 or len(b) == 0:
            continue
        statistic, p = _run_test(a, b, test)
        rows.append({
            "timepoint": timepoint,
            "time": group["time"].mean(),
            "green_intensity": intensity,
            "n_a": len(a), "n_b": len(b),
            "mean_a": np.mean(a), "mean_b": np.mean(b),
            "sd_a": np.std(a, ddof=1) if len(a) > 1 else np.nan,
            "sd_b": np.std(b, ddof=1) if len(b) > 1 else np.nan,
            "ratio_a_over_b": np.mean(a) / np.mean(b) if np.mean(b) != 0 else np.nan,
            "statistic": statistic,
            "p": p,
        })

    stats = pd.DataFrame(rows)
    if stats.empty:
        return stats
    stats["p_adj"] = _adjust_pvalues(stats["p"].to_numpy(), correction)
    stats["stars"] = [_p_to_stars(p) for p in stats["p_adj"]]
    stats["test"] = f"{test} ({group_a} vs {group_b}, correction={correction})"
    return stats


def _import_statsmodels():
    """Import statsmodels lazily (only needed for ANOVA); clear message if missing."""
    try:
        import statsmodels.api as sm
        import statsmodels.formula.api as smf
    except ImportError as error:
        raise ImportError("ANOVA needs statsmodels: pip install statsmodels") from error
    return sm, smf


def _fmt_p(p):
    """Format a p-value for display: '<0.001' or 2 significant figures."""
    if p is None or np.isnan(p):
        return "n/a"
    return "<0.001" if p < 0.001 else f"{p:.2g}"


def anova_stats(raw_dataframe, measurement, group_a="CS", group_b="KS",
                across="intensity", timepoints=None, intensities=None,
                log=True, correction="holm", replicate="auto",
                config=DefaultConfig):
    """
    Two-way ANOVA of ``group_a`` vs ``group_b`` with post-hoc comparisons.

    Only the two compared cell types enter the model (WT / media are very
    different and would distort the pooled variance).

    ``across="intensity"`` - one model per timepoint (panel):
        **two-way between-subjects ANOVA** ``value ~ strain * intensity``
        (Type II sums of squares). Valid because every well is an independent
        colony measured at one intensity.

    ``across="timepoint"`` - one model per intensity (panel):
        **mixed-design (repeated-measures) model** ``value ~ strain * time``
        with a random intercept per replicate, because the same wells are read
        at every timepoint. Terms are tested with Wald chi-squared tests.

    Post-hoc: ``group_a - group_b`` at every level of the x-axis factor,
    using the model's pooled error (more power than separate t-tests).
    Post-hoc p-values are corrected across ALL rows returned (the figure).

    :param raw_dataframe: long-format data (``sorted_data_df``)
    :param measurement: e.g. ``"GFP 395nm/OD600"``
    :param group_a: first cell type (default ``"CS"``)
    :param group_b: second cell type (default ``"KS"``)
    :param across: ``"intensity"`` or ``"timepoint"`` (see above)
    :param timepoints: int / list of read indices to include; None = all
    :param intensities: float / list of green intensities to include; None = all
    :param log: analyse log10(values) (default True). Recommended for
        fluorescence and ratio data: equalises variances and makes effects
        multiplicative, as the ANOVA assumes. Values must be > 0.
    :param correction: ``"holm"`` (default), ``"bh"`` or ``None``
    :param replicate: ``"auto"``, ``"well"`` or ``"colony"`` (see :func:`_replicate_values`)
    :param config: PlotConfig (exclusions)
    :return: ``(posthoc, omnibus)``:

        * ``posthoc`` - same columns as :func:`compare_conditions_stats`
          (``timepoint, time, green_intensity, n_a, n_b, mean_a, mean_b, ...,
          statistic, p, p_adj, stars, test``)
        * ``omnibus`` - one row per panel x term (strain, intensity/time,
          interaction) with ``statistic`` (F or chi2), ``df``, ``df_resid``, ``p``
    """
    sm, smf = _import_statsmodels()
    if across not in ("intensity", "timepoint"):
        raise ValueError("across must be 'intensity' or 'timepoint'")

    reps = _replicate_values(raw_dataframe, measurement, config, replicate)
    reps = reps[reps["cell_type"].isin([group_a, group_b])].copy()

    if timepoints is not None:
        timepoints = [timepoints] if np.isscalar(timepoints) else list(timepoints)
        reps = reps[reps["timepoint"].isin(timepoints)]
    if intensities is not None:
        intensities = [intensities] if np.isscalar(intensities) else list(intensities)
        reps = reps[np.isclose(reps["green_intensity"].to_numpy()[:, None], intensities).any(axis=1)]
    reps = reps.dropna(subset=["value"])

    if log:
        if (reps["value"] <= 0).any():
            raise ValueError("log=True needs all values > 0 (e.g. no negative blank-corrected values). "
                             "Use log=False.")
        reps["y"] = np.log10(reps["value"])
    else:
        reps["y"] = reps["value"]

    panel_col, x_col = (("timepoint", "green_intensity") if across == "intensity"
                        else ("green_intensity", "timepoint"))
    x_name = "intensity" if across == "intensity" else "time"

    omnibus_rows, post_rows = [], []
    for panel, d in reps.groupby(panel_col):
        d = d.copy()
        d["strain"] = d["cell_type"]
        d["x"] = d[x_col].map(str)                 # x-axis factor as categorical
        d["cell"] = d["strain"] + "|" + d["x"]     # one level per bar

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")        # convergence/boundary notes from MixedLM

            if across == "intensity":
                # --- Two-way between-subjects ANOVA ---------------------------
                model = smf.ols("y ~ C(strain) * C(x)", data=d).fit()
                table = sm.stats.anova_lm(model, typ=2)
                df_resid = table.loc["Residual", "df"]
                for term, label in [("C(strain)", "strain"), ("C(x)", x_name),
                                    ("C(strain):C(x)", f"strain x {x_name}")]:
                    omnibus_rows.append({panel_col: panel, "term": label, "test": "F",
                                         "statistic": table.loc[term, "F"],
                                         "df": table.loc[term, "df"], "df_resid": df_resid,
                                         "p": table.loc[term, "PR(>F)"]})
                # Cell-means model: same residuals, one coefficient per bar
                cell_model = smf.ols("y ~ C(cell) - 1", data=d).fit()
                names = list(cell_model.params.index)
            else:
                # --- Mixed model: strain between, time within (repeated) -----
                # Sum-to-zero coding so "strain" is the average effect over
                # time (a true main effect), not the effect at the first read
                model = smf.mixedlm("y ~ C(strain, Sum) * C(x, Sum)", d,
                                    groups=d["replicate_id"]).fit(reml=True)
                wald = model.wald_test_terms(scalar=True).table
                for term, label in [("C(strain, Sum)", "strain"), ("C(x, Sum)", x_name),
                                    ("C(strain, Sum):C(x, Sum)", f"strain x {x_name}")]:
                    omnibus_rows.append({panel_col: panel, "term": label, "test": "Wald chi2",
                                         "statistic": float(wald.loc[term, "statistic"]),
                                         "df": wald.loc[term, "df_constraint"], "df_resid": np.nan,
                                         "p": float(wald.loc[term, "pvalue"])})
                cell_model = smf.mixedlm("y ~ C(cell) - 1", d, groups=d["replicate_id"]).fit(reml=True)
                names = list(cell_model.fe_params.index)

        # --- Post-hoc contrasts: group_a - group_b at each x level -------------
        for xv in sorted(d[x_col].unique()):
            name_a = f"C(cell)[{group_a}|{xv}]"
            name_b = f"C(cell)[{group_b}|{xv}]"
            if name_a not in names or name_b not in names:
                continue
            contrast = np.zeros((1, len(names)))
            contrast[0, names.index(name_a)] = 1
            contrast[0, names.index(name_b)] = -1
            result = cell_model.t_test(contrast)

            a = d.loc[(d["strain"] == group_a) & (d[x_col] == xv), "value"]
            b = d.loc[(d["strain"] == group_b) & (d[x_col] == xv), "value"]
            row = {"timepoint": panel if panel_col == "timepoint" else xv,
                   "green_intensity": panel if panel_col == "green_intensity" else xv}
            row.update({
                "time": d.loc[d[x_col] == xv, "time"].mean() if panel_col == "green_intensity"
                        else d["time"].mean(),
                "n_a": len(a), "n_b": len(b),
                "mean_a": a.mean(), "mean_b": b.mean(),
                "sd_a": a.std(), "sd_b": b.std(),
                "ratio_a_over_b": a.mean() / b.mean() if b.mean() != 0 else np.nan,
                "statistic": float(np.squeeze(result.tvalue)),
                "p": float(np.squeeze(result.pvalue)),
            })
            post_rows.append(row)

    posthoc = pd.DataFrame(post_rows)
    omnibus = pd.DataFrame(omnibus_rows)
    if posthoc.empty:
        return posthoc, omnibus

    design = ("two-way ANOVA" if across == "intensity" else "mixed model (repeated measures)")
    posthoc["p_adj"] = _adjust_pvalues(posthoc["p"].to_numpy(), correction)
    posthoc["stars"] = [_p_to_stars(p) for p in posthoc["p_adj"]]
    posthoc["test"] = (f"{design}{' on log10' if log else ''} + post-hoc "
                       f"{group_a} vs {group_b} (correction={correction})")
    return posthoc, omnibus


def plot_comparison_bars(raw_dataframe,
                         measurement="GFP 395nm/OD600",
                         across="intensity",
                         timepoints=None,
                         intensities=None,
                         conditions=("CS", "KS", "WT", "media"),
                         normalise=None,
                         compare=("CS", "KS"),
                         test="welch",
                         anova_log=True,
                         stats_levels=None,
                         correction="holm",
                         p_format="stars",
                         replicate="auto",
                         label_map=None,
                         colors=None,
                         config=DefaultConfig,
                         save_image=False,
                         title_extra=""):
    """
    Grouped bar chart comparing cell types, across intensities or timepoints,
    with optional fold change and CS-vs-KS significance brackets.

    Layout:
        * ``across="intensity"`` - x-axis = green intensity; one panel per
          timepoint in ``timepoints`` (int or list).
        * ``across="timepoint"`` - x-axis = read (labelled in hours); one
          panel per intensity in ``intensities`` (float or list).
    Within each x group there is one bar per cell type in ``conditions``.

    Bars show mean +/- SD across replicates (see ``replicate``). With
    ``normalise="WT"`` or ``"media"`` each bar is divided by the reference
    mean at the same timepoint and intensity; the SD is propagated
    (relative errors in quadrature) and a dashed line marks 1. Where the
    reference is missing (e.g. media only at 0 and 1.4) the bar is left empty.

    Significance: ``compare=(a, b)`` adds a bracket over the two bars at
    every x, labelled with the corrected p-value from
    :func:`compare_conditions_stats` (family = everything in the figure).

    :param raw_dataframe: long-format data (``sorted_data_df``)
    :param measurement: e.g. ``"GFP 395nm/OD600"`` or ``"OD600"``
    :param across: ``"intensity"`` or ``"timepoint"``
    :param timepoints: across="intensity": panel timepoint(s) (0-based read
        index; use :func:`timepoint_at_hours` to convert from hours;
        None = last read).
        across="timepoint": optional subset of reads for the x-axis (None = all)
    :param intensities: across="timepoint": panel intensity/intensities
        (None = highest intensity). across="intensity": optional subset for
        the x-axis (None = all)
    :param conditions: cell types to draw, in bar order
    :param normalise: ``None`` (raw values), ``"WT"`` or ``"media"``
    :param compare: pair of cell types to test, or ``None`` for no stats
    :param test: ``"welch"`` (default), ``"student"``, ``"welch_log"``,
        ``"mannwhitney"`` - one test per bar pair; or ``"anova"`` - one model
        per panel (two-way ANOVA across intensities, repeated-measures mixed
        model across time) with post-hoc CS vs KS contrasts, see
        :func:`anova_stats`. ANOVA results are shown under each panel title.
    :param anova_log: ANOVA only - analyse log10 values (default True)
    :param stats_levels: optional subset of x-axis levels to include in the
        statistics (intensities if across="intensity", read indices if
        across="timepoint"). All bars are still drawn; only these get
        brackets and enter the ANOVA / correction. E.g. ``[0.28, 0.56, 1.4]``
        to test the responsive range only.
    :param correction: ``"holm"`` (default), ``"bh"`` or ``None``
    :param p_format: label above brackets - ``"stars"``, ``"p"`` or ``"both"``
    :param replicate: ``"auto"``, ``"well"`` or ``"colony"``
    :param label_map: optional ``{cell_type: readable name}`` for the legend,
        e.g. ``cell_label_map``
    :param colors: optional ``{cell_type: colour}``; default uses
        ``config.markercolor_map`` (media drawn grey)
    :param config: PlotConfig (exclusions, colours, save path)
    :param save_image: True -> save PNG and a stats CSV in ``figures/``
    :param title_extra: text appended to the title (and filename)
    :return: ``(summary, stats)`` - bar values per panel/x/cell type, and the
        stats table from :func:`compare_conditions_stats` (empty if compare=None).
        With ``test="anova"`` the omnibus table is in ``stats.attrs["anova"]``.
    """
    if across not in ("intensity", "timepoint"):
        raise ValueError("across must be 'intensity' or 'timepoint'")
    if normalise not in (None, "WT", "media"):
        raise ValueError("normalise must be None, 'WT' or 'media'")

    label_map = label_map or {}
    conditions = list(conditions)

    # --- 1. Replicate values and selection -----------------------------------
    reps = _replicate_values(raw_dataframe, measurement, config, replicate)

    def as_list(v):
        return None if v is None else ([v] if np.isscalar(v) else list(v))

    timepoints = as_list(timepoints)
    intensities = as_list(intensities)

    if across == "intensity":
        panel_col, x_col = "timepoint", "green_intensity"
        panels = timepoints if timepoints is not None else [int(reps["timepoint"].max())]
        if intensities is not None:
            reps = reps[np.isclose(reps["green_intensity"].to_numpy()[:, None], intensities).any(axis=1)]
        reps = reps[reps["timepoint"].isin(panels)]
    else:
        panel_col, x_col = "green_intensity", "timepoint"
        panels = intensities if intensities is not None else [float(reps["green_intensity"].max())]
        if timepoints is not None:
            reps = reps[reps["timepoint"].isin(timepoints)]
        reps = reps[np.isclose(reps["green_intensity"].to_numpy()[:, None], panels).any(axis=1)]

    if reps.empty:
        raise ValueError("No data for the selected timepoints / intensities.")

    # --- 2. Mean, SD, n per panel x x-value x cell type -----------------------
    summary = (reps.groupby(["timepoint", "green_intensity", "cell_type"], as_index=False)
                   .agg(time=("time", "mean"), mean=("value", "mean"),
                        std=("value", "std"), n=("value", "count")))

    # --- 3. Optional fold change vs WT or media -------------------------------
    if normalise is not None:
        ref = (summary[summary["cell_type"] == normalise]
               [["timepoint", "green_intensity", "mean", "std"]]
               .rename(columns={"mean": "ref_mean", "std": "ref_std"}))
        summary = summary.merge(ref, on=["timepoint", "green_intensity"], how="left")
        summary["ref_mean"] = summary["ref_mean"].replace(0, np.nan)
        fold = summary["mean"] / summary["ref_mean"]
        rel_err = np.sqrt((summary["std"].fillna(0) / summary["mean"]) ** 2
                          + (summary["ref_std"].fillna(0) / summary["ref_mean"]) ** 2)
        summary["bar"] = fold
        summary["bar_err"] = fold.abs() * rel_err
        ylabel = f"{measurement}\nfold change vs {normalise}"
    else:
        summary["bar"] = summary["mean"]
        summary["bar_err"] = summary["std"].fillna(0)
        ylabel = measurement

    # --- 4. Statistics ---------------------------------------------------------
    stats, omnibus = pd.DataFrame(), pd.DataFrame()
    stats_timepoints = sorted(reps["timepoint"].unique())
    stats_intensities = sorted(reps["green_intensity"].unique())
    if stats_levels is not None:
        if across == "intensity":
            stats_intensities = as_list(stats_levels)
        else:
            stats_timepoints = as_list(stats_levels)

    if compare is not None and test == "anova":
        # One model per panel; post-hoc p-values corrected across the figure
        stats, omnibus = anova_stats(
            raw_dataframe, measurement, compare[0], compare[1], across=across,
            timepoints=stats_timepoints,
            intensities=stats_intensities,
            log=anova_log, correction=correction, replicate=replicate, config=config)
        stats.attrs["anova"] = omnibus
    elif compare is not None:
        stats = compare_conditions_stats(
            raw_dataframe, measurement, compare[0], compare[1],
            timepoints=stats_timepoints,
            intensities=stats_intensities,
            test=test, correction=correction, replicate=replicate, config=config)

    # --- 5. Figure -------------------------------------------------------------
    # Colours: markercolor_map where defined, otherwise a distinct colour per
    # cell type (so missing entries are never all drawn in the same default)
    fallback = plt.get_cmap("tab10").colors
    default_colors = {c: config.markercolor_map.get(c, fallback[i % len(fallback)])
                      for i, c in enumerate(conditions)}
    if "media" in conditions:
        default_colors["media"] = "grey"   # keep media distinct from WT (both black in markercolor_map)
    colors = {**default_colors, **(colors or {})}

    n_bars = len(conditions)
    width = 0.8 / n_bars
    time_labels = summary.groupby("timepoint")["time"].mean()

    fig, axes = plt.subplots(len(panels), 1, squeeze=False, sharex=False,
                             figsize=(max(7, 1.3 * summary[x_col].nunique() * n_bars / 2.5 + 3),
                                      3.6 * len(panels)))
    axes = axes[:, 0]

    for ax, panel in zip(axes, panels):
        if panel_col == "green_intensity":
            data = summary[np.isclose(summary["green_intensity"], panel)]
            panel_title = f"Green intensity {panel:g}"
        else:
            data = summary[summary["timepoint"] == panel]
            panel_title = f"Read {panel} (~{time_labels.get(panel, np.nan):.1f} h)"

        # Add the ANOVA result for this panel as a second title line
        if not omnibus.empty:
            if panel_col == "green_intensity":
                terms = omnibus[np.isclose(omnibus["green_intensity"], panel)]
            else:
                terms = omnibus[omnibus["timepoint"] == panel]
            stat_name = "F" if (terms["test"] == "F").all() else "χ²"
            parts = [f"{r.term}: {stat_name}={r.statistic:.1f}, p{'' if _fmt_p(r.p).startswith('<') else '='}{_fmt_p(r.p)}"
                     for r in terms.itertuples()]
            panel_title += "\n" + ("ANOVA" if stat_name == "F" else "Mixed model") \
                           + " - " + ";  ".join(parts)

        x_values = sorted(data[x_col].unique())
        x = np.arange(len(x_values))
        bar_tops = {}   # (x index, cell type) -> top of bar + error, for brackets

        for i, cell in enumerate(conditions):
            subset = data[data["cell_type"] == cell].set_index(x_col).reindex(x_values)
            y = subset["bar"].to_numpy(dtype=float)
            err = subset["bar_err"].to_numpy(dtype=float)
            positions = x + (i - (n_bars - 1) / 2) * width
            ax.bar(positions, y, width, yerr=err, capsize=3,
                   color=colors.get(cell, "grey"), alpha=0.8,
                   edgecolor="black", linewidth=0.5,
                   label=label_map.get(cell, cell))
            for xi, (pos, top) in enumerate(zip(positions, np.nan_to_num(y) + np.nan_to_num(err))):
                bar_tops[(xi, cell)] = (pos, top)

        # Significance brackets between the two compared cell types
        if compare is not None and not stats.empty and all(c in conditions for c in compare):
            y_max = max((t for _, t in bar_tops.values()), default=1) or 1
            step = 0.05 * y_max
            highest = y_max
            for xi, xv in enumerate(x_values):
                if panel_col == "green_intensity":
                    row = stats[(stats["timepoint"] == xv) & np.isclose(stats["green_intensity"], panel)]
                else:
                    row = stats[(stats["timepoint"] == panel) & np.isclose(stats["green_intensity"], xv)]
                if row.empty or np.isnan(row["p_adj"].iloc[0]):
                    continue
                p_adj = row["p_adj"].iloc[0]
                stars = row["stars"].iloc[0]
                if p_format == "p":
                    text = f"p={p_adj:.2g}"
                elif p_format == "both":
                    text = f"{stars}\np={p_adj:.2g}"
                else:
                    text = stars

                (xa, ta), (xb, tb) = bar_tops[(xi, compare[0])], bar_tops[(xi, compare[1])]
                y0 = max(ta, tb) + step
                ax.plot([xa, xa, xb, xb], [y0, y0 + step / 2, y0 + step / 2, y0],
                        color="black", linewidth=0.8)
                ax.text((xa + xb) / 2, y0 + step * 0.6, text,
                        ha="center", va="bottom", fontsize=8)
                highest = max(highest, y0 + step * (3 if p_format == "both" else 2))
            ax.set_ylim(top=highest * 1.05)

        # Axis formatting
        ax.set_xticks(x)
        if x_col == "green_intensity":
            ax.set_xticklabels([f"{v:g}" for v in x_values])
            ax.set_xlabel("Green light intensity")
        else:
            ax.set_xticklabels([f"{time_labels.get(v, np.nan):.0f} h" for v in x_values])
            ax.set_xlabel("Time")
        ax.set_ylabel(ylabel)
        ax.set_title(panel_title, fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.2, linestyle="--")
        if normalise is not None:
            ax.axhline(1, color="grey", linestyle="--", linewidth=1)

    # Legend once, to the right of the top panel
    axes[0].legend(title="Cell type", loc="upper left", bbox_to_anchor=(1.01, 1.0),
                   borderaxespad=0.0, fontsize="small", title_fontsize="small")

    # Title and footnote describing the test
    norm_text = f" - fold vs {normalise}" if normalise else ""
    title = f"{measurement} bars across {across}{norm_text} - {title_extra}"
    fig.suptitle(title, fontsize=11)
    if compare is not None and not stats.empty:
        n_text = f"n = {int(stats['n_a'].min())}-{int(stats['n_a'].max())} ({compare[0]}), " \
                 f"{int(stats['n_b'].min())}-{int(stats['n_b'].max())} ({compare[1]}) replicates"
        if test == "anova":
            design = ("Two-way ANOVA (strain x intensity, Type II) per panel"
                      if across == "intensity" else
                      "Mixed model (strain x time, random intercept per well, Wald tests) per panel")
            method = (f"{design}{' on log10 values' if anova_log else ''}, {compare[0]} and "
                      f"{compare[1]} only. Brackets: post-hoc {compare[0]} vs {compare[1]} contrasts "
                      f"from the model, {correction or 'no'} correction across {len(stats)} comparisons")
        else:
            method = (f"{compare[0]} vs {compare[1]}: {test} test, {correction or 'no'} correction "
                      f"across {len(stats)} comparisons")
        fig.text(0.01, -0.01,
                 f"{method}. {n_text}. "
                 "* p<0.05, ** p<0.01, *** p<0.001, ns not significant. Bars: mean ± SD.",
                 fontsize=7, ha="left", va="top", wrap=True)
    fig.tight_layout()

    # --- 6. Save or show (+ stats table) ---------------------------------------
    if save_image and not stats.empty:
        try:
            folder = Path(config.save_file_path) / "figures"
            folder.mkdir(parents=True, exist_ok=True)
            stats.to_csv(folder / f"{title.replace('/', '_div_')}_stats.csv", index=False)
            if not omnibus.empty:
                omnibus.to_csv(folder / f"{title.replace('/', '_div_')}_anova.csv", index=False)
        except OSError as error:
            print(f"Could not save stats table: {error}")
    save_or_show(fig, title, config, save_image)

    return summary, stats


if __name__ == "__main__":
    main(cell_label_map, markerstyle_map, markercolor_map, linestyle_map, line_color_map, plot_exclude)
