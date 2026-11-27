"""
legend_helpers.py
=================

Drop-in helpers for plate_reader_analysis_*.py that

1. remove excluded conditions from a dataframe in one place
   (:func:`apply_exclusions`), and
2. draw legends that list ONLY the values actually plotted, stacked outside
   the axes on the right so they never overlap each other or the plot
   (:func:`add_legends`).

Usage inside a plot function (replaces the three ``axs.legend(...)`` blocks)::

    raw_data     = apply_exclusions(raw_data, config)
    summary_data = apply_exclusions(summary_data, config)
    ... plotting ...
    axs.set_xlabel(...); axs.set_ylabel(...); axs.set_title(...)
    fig.tight_layout()                         # call BEFORE add_legends
    add_legends(axs, summary_data, config)     # last thing before saving

and save with ``bbox_inches="tight"`` so nothing is clipped::

    fig.savefig(path, bbox_inches="tight")
"""

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


# Default columns, legend titles and order of the legends (top to bottom)
LEGEND_SPECS = {
    "cells": "Cell type",
    "media": "Media",
    "green_intensity": "Green light intensity",
}


def apply_exclusions(df, config, columns=("cells", "media", "green_intensity", "red_intensity")):
    """
    Return a copy of ``df`` without the rows listed in ``config.plot_exclude``.

    Uses the same rule as the existing ``if ... in config.plot_exclude[...]``
    checks (exact match), so plotted data and legends always agree.

    :param df: long-format or summary dataframe
    :param config: PlotConfig with a ``plot_exclude`` dict
    :param columns: which exclusion keys to apply. ``plot_by_intensity`` only
        applies ``("cells", "media")`` because intensity is its x-axis.
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
    Build legend handles for one column.

    Missing style entries fall back to the ``default_*`` values in PlotConfig
    instead of raising KeyError.
    """
    label_map = label_map or {}
    handles = []
    for value in values:
        label = label_map.get(value, f"{value:g}" if isinstance(value, (int, float, np.number)) else str(value))

        if column == "cells":
            colour = config.markercolor_map.get(value, config.default_markercolor)
            handles.append(Line2D([], [], linestyle="None",
                                  marker=config.markerstyle_map.get(value, config.default_markerstyle),
                                  markerfacecolor=colour, markeredgecolor=colour,
                                  markersize=6, label=label))
        elif column == "media":
            handles.append(Line2D([], [], color="black", linewidth=2,
                                  linestyle=config.linestyle_map.get(value, config.default_linestyle),
                                  label=label))
        elif column == "green_intensity":
            handles.append(Line2D([], [], linewidth=2,
                                  color=config.line_color_map.get(value, config.default_linecolor),
                                  label=label))
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
    resizing), and save with ``bbox_inches="tight"``.

    :param ax: matplotlib Axes that was plotted on
    :param data: the dataframe that was actually plotted (already filtered with
        :func:`apply_exclusions`)
    :param config: PlotConfig with the style maps
    :param legends: dict ``{column: title}`` in display order. Default:
        cells, media, green intensity. E.g. for ``plot_by_intensity`` pass
        ``{"cells": "Cell type", "media": "Media"}``.
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
