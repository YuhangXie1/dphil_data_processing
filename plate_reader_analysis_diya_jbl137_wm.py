import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import copy
import json
import os
from pathlib import Path
from datetime import datetime
from os import listdir, path
import re
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from enum import Enum
from dataclasses import dataclass
from typing import Dict, Any, Optional

#Defining default config
@dataclass
class PlotConfig:
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
    return PlotConfig(
        markerstyle_map = {"JBL001":"s", "YX002":"*", "YX001":"o", "media":"^"},
        markercolor_map = {"JBL001":"black", "YX002":"orange", "YX001": "blue", "media":"black"},
        linestyle_map = {"WM-met+": "solid", "WM-met-": "dashed", "LB": "solid"},
        line_color_map = {},  # build from data later
        line_alpha_map = {},
        plot_exclude = {"cells":[], "media":[], "green_intensity":[], "red_intensity":[]},
        alpha_used = False,
        save_file_path = save_file_path
    )
DefaultConfig = default_plot_config()

### FUNCTIONS ###

def mean_std_cv(data):
    #data should be array containing 1 or more arrays
    mean = np.mean(data, axis = 0)
    std_dev = np.std(data, axis = 0)
    cv = (std_dev/mean)*100
    return [mean, std_dev, cv]

def process_data_to_mean(data_dict):
    data_dict_repeats = {}
    for key, value in data_dict.items():
        #if key == "media":
        #    data_dict_repeats["media"].append(value)
        #parsing name
        namemet_repeat = key.rsplit("_",1)
        if namemet_repeat[0] in data_dict_repeats.keys():
            data_dict_repeats[namemet_repeat[0]].append(value)
        else:
            data_dict_repeats[namemet_repeat[0]] = [value]

    print(data_dict_repeats)
    data_dict_processed = {}
    for key, value in data_dict_repeats.items():
        #getting mean, std, cv
        mean_std_cv_array = mean_std_cv(value) #mean_std_cv_array = [[mean_od,mean_gfp, mean_gfp/od],[std_od,std_gfp, std_gfp/od],[cv_od,cv_gfp,cv_gfp/od]]

        #new dictionary with processed data. value = [mean, std, cv]
        data_dict_processed[key] = mean_std_cv_array
    print(data_dict_processed)
    return data_dict_processed

def build_intensity_colour_map(green_intensities, cmap_name = None) -> Dict[Any,tuple]:
    green_intensities_set = set(green_intensities)
    min_intensity = min(green_intensities_set)
    max_intensity = max(green_intensities_set)

    if cmap_name is None:
        cmap = LinearSegmentedColormap.from_list("red_to_green",[(1,0,0), (0,1,0)])
    else:
        cmap = plt.get_cmap(cmap_name)

    if max_intensity == min_intensity:
        print("only one green intensity")
        return {}
    
    norm = mcolors.Normalize(vmin=min_intensity,vmax=max_intensity)

    return {
            intensity: cmap(norm(intensity))
            for intensity in green_intensities_set
        }

def create_plot_handles(config: PlotConfig):
    #legend handles
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

def load_data(filepath):
    
    #initialises columns of dataframe
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

    #Loading data
    data = pd.read_csv(filepath, header = [0,1], index_col= 0)
   
    #Rebuild the MultiIndex with proper Timestamp objects
    timestamps = pd.to_datetime(data.columns.get_level_values("timestamp"), dayfirst=True)
    measurements = data.columns.get_level_values("measurement")
    data.columns = pd.MultiIndex.from_arrays([measurements, timestamps], names=["measurement", "timestamp"])
    rows = []

    #Populating the DF
    for well in data.index:
        
        #Plate settings
        cells = plate_map[well][0]
        media = plate_map[well][1]
        green_intensity = plate_map[well][2]
        red_intensity = plate_map[well][3]

        #Populate for each mode and timestamp
        for measurement, timestamp in data.columns:
            timestamps = data[measurement].columns
            timepoint = timestamps.get_loc(timestamp)
            initial_time = timestamps.min()
            time = (timestamp - initial_time).total_seconds() /3600
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

    #Calculating signal div OD600
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

                #finding corresponding OD row
                od_row = select_od_data[(select_od_data["well"] == well) & (select_od_data["timepoint"] == timepoint)]
                od_value = od_row["value"].iloc[0]

                #calculating ratio
                measure_value = measure_row["value"]
                ratio = measure_value / od_value

                #writing new data row
                new_row = measure_row.copy()
                new_row["measurement"] = calculated_measurement
                new_row["value"] = ratio
                new_rows.append(new_row)

    calculated_data_df = pd.DataFrame(new_rows)
    sorted_data_df = pd.concat([sorted_data_df, calculated_data_df], ignore_index= True)

    #finding means and averages
    summary_df = sorted_data_df.groupby(["cells","media","green_intensity","red_intensity","timestamp","time","timepoint","measurement"], as_index=False).agg(
        min_timestamp = ("timestamp", "min"),
        mean = ("value","mean"),
        std = ("value","std"),
        count = ("value","count")
    )

    # Generating color map
    DefaultConfig.line_color_map = build_intensity_colour_map(sorted_data_df["green_intensity"].unique())

    return sorted_data_df, summary_df

def plot_timecourse(raw_dataframe, summary_dataframe, measurement, plot_type, ylabel: str | None = None, title: str | None = None, title_extra: str = "",
                            xlabel = "Time (hrs)",
                            config: PlotConfig = DefaultConfig, save_image = False,
                           ):
    """
    Plots y_data over time
    
    :param dataframe: dataframe where the data is (raw)
    :param summary_dataframe: dataframe where the average, std, count data is
    :param measurement: data type to plot (e.g. OD600 or GFP 395nm)
    :param plot_type: select whether 'average' are plotted or 'all' data
    :param ylabel: optional override for the ylabel text
    :type ylabel: str | None
    :param title: optional override for title text
    :type title: str | None
    :param title_extra: optional extra to tag onto the end of the title text
    :type title_extra: str
    :param xlabel: optional override for the xlabel text
    :param config: config file for styling and saving the image
    :type config: PlotConfig
    :param save_image: if True, then image is saved instead of displayed
    """

    if ylabel is None:
        ylabel = measurement

    if title is None:
        title = f"{ylabel} - {plot_type} - {title_extra}"
    else:
        title = title + " - " + title_extra

    #slicing data
    summary_data = summary_dataframe.loc[summary_dataframe["measurement"] == measurement].copy() 
    raw_data = raw_dataframe.loc[raw_dataframe["measurement"] == measurement].copy() 
    

    fig, axs = plt.subplots()
    if plot_type == "all" or plot_type == "both":
        if plot_type == "both":
            select_alpha = 0.2
        else:
            select_alpha = 1.0

        for well, group in raw_data.groupby(["well"]):
            row = group.iloc[0]
            cells = row["cells"]
            media = row["media"]
            green_intensity = row["green_intensity"]
            red_intensity = row["red_intensity"]

            if (cells in config.plot_exclude["cells"]
                        or media in config.plot_exclude["media"]
                        or green_intensity in config.plot_exclude["green_intensity"]
                        or red_intensity in config.plot_exclude["red_intensity"]):
                    continue
            
            axs.plot(group["time"], group["value"],

                        color = config.line_color_map[green_intensity],
                        marker = config.markerstyle_map[cells],
                        markerfacecolor = config.markercolor_map[cells],
                        markeredgecolor = config.markercolor_map[cells], markersize = 3.0,
                        linestyle = config.linestyle_map[media], linewidth = 1.0,
                        alpha = select_alpha,
                        )  

    if plot_type == "average" or plot_type == "both":
        for (cells, media, green_intensity, red_intensity), group in summary_data.groupby(["cells","media","green_intensity","red_intensity"]):
            if (cells in config.plot_exclude["cells"]
                        or media in config.plot_exclude["media"]
                        or green_intensity in config.plot_exclude["green_intensity"]
                        or red_intensity in config.plot_exclude["red_intensity"]):
                    continue

            axs.errorbar(group["time"], group["mean"],
                        yerr = group["std"], capsize = 2.0,

                        color = config.line_color_map[green_intensity],
                        marker = config.markerstyle_map[cells],
                        markerfacecolor = config.markercolor_map[cells],
                        markeredgecolor = config.markercolor_map[cells], markersize = 3.0,
                        linestyle = config.linestyle_map[media], linewidth = 1.0,
                        alpha = 1.0,
                        )

    cell_handles, intensity_handles, media_handles = create_plot_handles(config)

    leg1 = axs.legend(
        handles=cell_handles,
        title="Cell type",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00),
        borderaxespad=0.0,
    )
    axs.add_artist(leg1)

    leg2 = axs.legend(
        handles=media_handles,
        title="Media",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.65),
        borderaxespad=0.0,
    )
    axs.add_artist(leg2)

    leg3 = axs.legend(
        handles=intensity_handles,
        title="Green light intensity",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.35),
        borderaxespad=0.0,
    )

    axs.add_artist(leg3)
    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_title(title)
    fig.tight_layout(rect=[0, 0, 0.75, 1])

    if save_image is True:
        try:
            Path(os.path.join(config.save_file_path, "figures")).mkdir(parents = True, exist_ok = True)
            save_title = title.replace("/","_div_")
            plt.savefig(os.path.join(config.save_file_path, "figures", f"{save_title}.png"))
            plt.close(fig)
        except:
            print("Could not save image, filepath not valid")
    else:
        plt.show()
        plt.close(fig)



def plot_by_intensity(raw_dataframe, summary_dataframe, measurement, plot_type, timepoints: int | list[int], ylabel: str | None = None, title: str | None = None,
                              title_extra: str = "",
                              xlabel = "Green light intensity %",
                                config: PlotConfig = DefaultConfig, save_image = False):
    """
    Plots y_data over intensity at a chosen time point
    
    :param raw_dataframe: dataframe where the raw data is
    :param summary_dataframe: dataframe where the summary data is
    :param measurement: the measurement to plot (e.g. OD600 or GFP 395nm)
    :param plot_type: select whether 'average' are plotted or 'all' data
    :param timepoint: chosen time point (number of reading, e.g. 0 = first reading, 1 = second reading, etc.) or series of timepoints
    :param ylabel: optional override for the ylabel text
    :type ylabel: str | None
    :param title: optional override for title text
    :type title: str | None
    :param title_extra: optional extra to tag onto the end of the title text
    :type title_extra: str
    :param xlabel: optional override for the xlabel text
    :param config: config file for styling and saving the image
    :type config: PlotConfig
    :param save_image: if True, then image is saved instead of displayed
    """
    if ylabel is None:
        ylabel = measurement

    if title is None:
        title = f"{ylabel} - by intensity - {plot_type} - {title_extra}"
    else:
        title = title + " - " + title_extra

    if isinstance(timepoints, int): #single timepoint allows gradient colours to be used
        timepoints = [timepoints]
        colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
        cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)
        timepoints_is_array = False
    else:
        cmap = plt.get_cmap("viridis")  # Use a different colormap for multiple timepoints
        timepoints_is_array = True
        #TODO - implement different line colours for multiple timepoints

    fig, axs = plt.subplots()
    for timepoint in timepoints:
        #slicing data
        raw_data = raw_dataframe.loc[(raw_dataframe["measurement"] == measurement) & (raw_dataframe["timepoint"] == timepoint)].copy()
        summary_data = summary_dataframe.loc[(summary_dataframe["measurement"] == measurement) & (summary_dataframe["timepoint"] == timepoint)].copy()

        raw_data["green_intensity_percentage"] = raw_data["green_intensity"].apply(lambda x: 100*x/2.8)
        raw_data["red_intensity_percentage"] = raw_data["red_intensity"].apply(lambda x: 100*x/2.8)
        summary_data["green_intensity_percentage"] = summary_data["green_intensity"].apply(lambda x: 100*x/2.8)
        summary_data["red_intensity_percentage"] = summary_data["red_intensity"].apply(lambda x: 100*x/2.8)

        select_alpha = 1.0
        if plot_type == "all" or plot_type == "both":
            if plot_type == "both":
                select_alpha = 0.2

            for well, group in raw_data.groupby(["well"]):
                row = group.iloc[0]
                cells = row["cells"]
                media = row["media"]
                green_intensity = row["green_intensity"]
                red_intensity = row["red_intensity"]

                if (cells in config.plot_exclude["cells"]
                        or media in config.plot_exclude["media"]
                        or green_intensity in config.plot_exclude["green_intensity"]
                        or red_intensity in config.plot_exclude["red_intensity"]):
                    continue

                x = group["green_intensity_percentage"].values.copy()
                y = group["value"].values.copy()

                # move final point (green=2.8 & red=0) to the right
                is_final = ((group["green_intensity"] == 2.8) &(group["red_intensity"] == 0))
                x[is_final.values] = x.max() + 10  # move to the right
                
                order = np.argsort(x)
                x = x[order]
                y = y[order]

                # No line segments as points are not connected to each other
                norm = plt.Normalize(0, 100)
                # Add scatter points with gradient colors
                scatter = axs.scatter(x, y, c=row["green_intensity_percentage"], cmap=cmap, norm=norm, s=20, zorder=5,
                                    marker = config.markerstyle_map[cells],
                                    facecolor = config.markercolor_map[cells],
                                    edgecolor = config.markercolor_map[cells],
                                    alpha = select_alpha,
                                    )
                
        if plot_type == "average" or plot_type == "both":
                for (cells, media), group in summary_data.groupby(["cells","media"]):
                    if (cells in config.plot_exclude["cells"]
                                or media in config.plot_exclude["media"]):
                            continue
        
                    x = group["green_intensity_percentage"].values.copy()
                    y = group["mean"].values.copy()
                    yerr = group["std"].values.copy()
        
                    # move final point (green=2.8 & red=0) to the right
                    is_final = ((group["green_intensity"] == 2.8) &(group["red_intensity"] == 0))
                    x[is_final.values] = x.max() + 10  # move to the right
                    
                    order = np.argsort(x)
                    x = x[order]
                    y = y[order]
                    yerr = yerr[order]
        
                    # Create line segments for gradient
                    points = np.array([x, y]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]], axis=1)
        
                    # Normalize x values for colormap
                    norm = plt.Normalize(x.min(), x.max())
                    lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
                    lc.set_linestyle(config.linestyle_map[media])
                    lc.set_alpha(select_alpha)
                    lc.set_array(x)
                    axs.add_collection(lc)
        
                    # Add scatter points with gradient colors
                    scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5,
                                        marker = config.markerstyle_map[cells],
                                        facecolor = config.markercolor_map[cells],
                                        edgecolor = config.markercolor_map[cells],
                                        alpha = select_alpha,
                                        )
                    
                    axs.errorbar(x, y, yerr=yerr, fmt='none',
                                ecolor= config.markercolor_map[cells], alpha=1.0, capsize=2.0,
                                )
    
    # labels
    cell_handles, intensity_handles, media_handles = create_plot_handles(config)
    leg1 = axs.legend(
        handles=cell_handles,
        title="Cell type",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00),
        borderaxespad=0.0,
    )
    axs.add_artist(leg1)

    leg2 = axs.legend(
        handles=media_handles,
        title="Media",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.65),
        borderaxespad=0.0,
    )
    axs.add_artist(leg2)

    x_axis = axs.set_xlabel(xlabel)
    x_axis.set_color("green")

    axs.set_ylabel(ylabel)
    axs.set_title(title)

    fig.tight_layout(rect=[0, 0, 0.75, 1])
    if save_image == True:
        try:
            Path(os.path.join(config.save_file_path, "figures")).mkdir(parents = True, exist_ok = True)
            save_title = title.replace("/","_div_")
            plt.savefig(os.path.join(config.save_file_path, "figures", f"{save_title}.png"))
            plt.close(fig)
        except:
            print("Could not save image, filepath not valid")
    else:
        plt.show()
        plt.close(fig)

def plot_by_intensity_foldchange(summary_dataframe, measurement, timepoints: int | list[int], ylabel: str | None = None, title: str | None = None,
                              title_extra: str = "",
                              xlabel = "Green light intensity %",
                                config: PlotConfig = DefaultConfig, save_image = False):
    """
    Plots fluorescence foldchange over intensity at a chosen time point. Fold change is calculated by the select green light divided by the 100% red light value
    
    :param summary_dataframe: dataframe where the summary data is
    :param measurement: the measurement to plot (e.g. OD600 or GFP 395nm)
    :param plot_type: select whether 'average' are plotted or 'all' data
    :param timepoints: chosen time point (number of reading, e.g. 0 = first reading, 1 = second reading, etc.) or series of timepoints
    :param ylabel: optional override for the ylabel text
    :type ylabel: str | None
    :param title: optional override for title text
    :type title: str | None
    :param title_extra: optional extra to tag onto the end of the title text
    :type title_extra: str
    :param xlabel: optional override for the xlabel text
    :param config: config file for styling and saving the image
    :type config: PlotConfig
    :param save_image: if True, then image is saved instead of displayed
    """
    if ylabel is None:
        ylabel = measurement

    if title is None:
        title = f"{ylabel} - by intensity - average - {title_extra}"
    else:
        title = title + " - " + title_extra

    if isinstance(timepoints, int): #single timepoint allows gradient colours to be used
        timepoints = [timepoints]
        colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
        cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)
        timepoints_is_array = False
    else:
        cmap = plt.get_cmap("viridis")  # Use a different colormap for multiple timepoints
        timepoints_is_array = True
        #TODO - implement different line colours for multiple timepoints

    fig, axs = plt.subplots()
    for timepoint in timepoints:
        #slicing data
        summary_data = summary_dataframe.loc[(summary_dataframe["measurement"] == measurement) & (summary_dataframe["timepoint"] == timepoint)].copy()
        summary_data["green_intensity_percentage"] = summary_data["green_intensity"].apply(lambda x: 100*x/2.8)
        summary_data["red_intensity_percentage"] = summary_data["red_intensity"].apply(lambda x: 100*x/2.8)

        for (cells, media), group in summary_data.groupby(["cells","media"]):
            if (cells in config.plot_exclude["cells"]
                        or media in config.plot_exclude["media"]):
                    continue

            green_data = group.loc[group["green_intensity"] > 0].copy()
            x = green_data["green_intensity_percentage"].values

            y_green = green_data["mean"].values
            yerr_green = green_data["std"].values
            y_red = group.loc[(group["red_intensity"] == 2.8) & (group["green_intensity"] == 0), "mean"].values.copy()
            yerr_red = group.loc[(group["red_intensity"] == 2.8) & (group["green_intensity"] == 0), "std"].values.copy()

            y = y_green / y_red  # Calculate fold change
            yerr = np.sqrt((yerr_green / y_green)**2 + (yerr_red / y_red)**2) * y  # Propagate error

            # move final point (green=2.8 & red=0) to the right
            is_final = ((green_data["green_intensity"] == 2.8) &(green_data["red_intensity"] == 0))
            x[is_final.values] = x.max() + 10  # move to the right
            
            order = np.argsort(x)
            x = x[order]
            y = y[order]
            yerr = yerr[order]

            # Create line segments for gradient
            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            # Normalize x values for colormap
            norm = plt.Normalize(x.min(), x.max())
            lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0, alpha=1.0)
            lc.set_linestyle(config.linestyle_map[media])
            lc.set_array(x)
            axs.add_collection(lc)

            # Add scatter points with gradient colors
            scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5,
                                marker = config.markerstyle_map[cells],
                                facecolor = config.markercolor_map[cells],
                                edgecolor = config.markercolor_map[cells],
                                alpha = 1.0,
                                )
            
            axs.errorbar(x, y, yerr=yerr, fmt='none',
                        ecolor= config.markercolor_map[cells], alpha=1.0, capsize=2.0,
                        )
    
    # labels
    cell_handles, intensity_handles, media_handles = create_plot_handles(config)
    leg1 = axs.legend(
        handles=cell_handles,
        title="Cell type",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00),
        borderaxespad=0.0,
    )
    axs.add_artist(leg1)

    leg2 = axs.legend(
        handles=media_handles,
        title="Media",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.65),
        borderaxespad=0.0,
    )
    axs.add_artist(leg2)

    x_axis = axs.set_xlabel(xlabel)
    x_axis.set_color("green")

    axs.set_ylabel(ylabel)
    axs.set_title(title)

    fig.tight_layout(rect=[0, 0, 0.75, 1])
    if save_image == True:
        try:
            Path(os.path.join(config.save_file_path, "figures")).mkdir(parents = True, exist_ok = True)
            save_title = title.replace("/","_div_")
            plt.savefig(os.path.join(config.save_file_path, "figures", f"{save_title}.png"))
            plt.close(fig)
        except:
            print("Could not save image, filepath not valid")
    else:
        plt.show()
        plt.close(fig)

def plot_timecourse_custom(dataframe, y_data, plot_type, ylabel: str | None = None, title: str | None = None, title_extra: str = "",
                            xlabel = "Time (hrs)",
                            config: PlotConfig = DefaultConfig, save_image = False,
                            row_filter = None,
                           ):
    """
    Plots y_data over time
    
    :param dataframe: dataframe where the data is
    :param y_data: data to plot
    :param plot_type: select whether 'average' are plotted or 'all' data
    :param ylabel: optional override for the ylabel text
    :type ylabel: str | None
    :param title: optional override for title text
    :type title: str | None
    :param title_extra: optional extra to tag onto the end of the title text
    :type title_extra: str
    :param xlabel: optional override for the xlabel text
    :param config: config file for styling and saving the image
    :type config: PlotConfig
    :param save_image: if True, then image is saved instead of displayed
    """

    if ylabel is None:
        ylabel = y_data

    if title is None:
        title = f"{ylabel} - {plot_type} - {title_extra}"
    else:
        title = title + " - " + title_extra


    fig, axs = plt.subplots()
    for row_index, row in dataframe.iterrows():

        # Bespoke filter goes here
        if row_filter is not None and not row_filter(row):
            continue

        index_name_cells = row["cells"]
        index_name_media = row["media"]
        index_name_green_intensity = row["green_intensity"]
        index_name_red_intensity = row["red_intensity"]

        if (index_name_cells in config.plot_exclude["cells"]
            or index_name_media in config.plot_exclude["media"]
            or index_name_green_intensity in config.plot_exclude["green_intensity"]
            or index_name_red_intensity in config.plot_exclude["red_intensity"]):
            continue

        if config.alpha_used is True:
            alpha = config.line_alpha_map[index_name_green_intensity]
        else:
            alpha = 1

        if plot_type == "average":
            axs.errorbar(row[f"{y_data}_timepoints"], row[f"{y_data}_average"],
                        yerr = row[f"{y_data}_std"], capsize = 2.0,

                        color = config.line_color_map[index_name_green_intensity],
                        marker = config.markerstyle_map[index_name_cells],
                        markerfacecolor = config.markercolor_map[index_name_cells],
                        markeredgecolor = config.markercolor_map[index_name_cells], markersize = 3.0,
                        linestyle = config.linestyle_map[index_name_media], linewidth = 1.0,
                        alpha = alpha,
                        )
            
        elif plot_type == "all":
            for repeat in row[f"{y_data}_raw_array"]:
                axs.plot(row[f"{y_data}_timepoints"], repeat,
                        
                            color = config.line_color_map[index_name_green_intensity],
                            marker = config.markerstyle_map[index_name_cells],
                            markerfacecolor = config.markercolor_map[index_name_cells],
                            markeredgecolor = config.markercolor_map[index_name_cells], markersize = 3.0,
                            linestyle = config.linestyle_map[index_name_media], linewidth = 1.0,
                            alpha = alpha,
                            )
        
        elif plot_type == "both":

            if row["cells"] == "JBL137":
                custom_line_color = {0.28: "blue", 0.0: "red", 2.8: "blue", 0.028:"blue"}
            else:
                custom_line_color = {0.28: "black", 0.0: "black", 2.8: "black", 0.028:"blue"}

            axs.errorbar(row[f"{y_data}_timepoints"], row[f"{y_data}_average"],
                        yerr = row[f"{y_data}_std"], capsize = 2.0,

                        color = custom_line_color[index_name_green_intensity],
                        marker = config.markerstyle_map[index_name_cells],
                        markerfacecolor = custom_line_color[index_name_green_intensity],
                        markeredgecolor = custom_line_color[index_name_green_intensity], markersize = 3.0,
                        linestyle = config.linestyle_map[index_name_media], linewidth = 1.0,
                        alpha = 1.0,
                        )
            

            for repeat in row[f"{y_data}_raw_array"]:
                if row["cells"] == "JBL001" or row["cells"] == "media":
                    continue
                axs.plot(row[f"{y_data}_timepoints"], repeat,
                        
                            color = custom_line_color[index_name_green_intensity],
                            linestyle = config.linestyle_map[index_name_media], linewidth = 1.0,
                            alpha = 0.2,
                            )
                            

    leg1 = axs.legend(
        handles=cell_handles,
        title="Cell type",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00),
        borderaxespad=0.0,
    )
    axs.add_artist(leg1)

    leg2 = axs.legend(
        handles=media_handles,
        title="Media",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.65),
        borderaxespad=0.0,
    )
    axs.add_artist(leg2)

    leg3 = axs.legend(
        handles=intensity_handles,
        title="Green light intensity",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.35),
        borderaxespad=0.0,
    )

    axs.add_artist(leg3)
    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_title(title)
    fig.tight_layout(rect=[0, 0, 0.75, 1])

    if save_image is True:
        try:
            Path(os.path.join(config.save_file_path, "figures")).mkdir(parents = True, exist_ok = True)
            save_title = title.replace("/","_div_")
            plt.savefig(os.path.join(config.save_file_path, "figures", f"{save_title}.png"))
            plt.close(fig)
        except:
            print("Could not save image, filepath not valid")
    else:
        plt.show()
        plt.close(fig)

def populate_plate_map(optical_power_map, cell_map, media_map):
    pass

key_rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
key_columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
key_wells = [str(row)+str(col) for row in key_rows for col in key_columns]
plate_map = {key : [] for key in key_wells}

NCOLS = 12
NROWS = 8
### MAIN ###

#copy and paste the light array

DEFAULT_OPTICAL_POWER = np.array([

	# Channel 0 (Color 0 or 4). Blue on v0.4c
	[[0 / (2**(7-row)) for col in range(NCOLS)] for row in range(NROWS)],

    # Channel 1 – Green light (checkerboard starting green at B2)
    np.array([
        # Columns 1–12 (A–L), Rows A–H
        [2.8,	0.56,   1.4,    2.8,	0.56,   1.4,    2.8,	0.56,   1.4,   	2.8,	0.56,   1.4,],  # Row A
        [0.0,  	0.28,  	0.28,  	0.0,  	0.28,  	0.28,	0.0,  	0.28,  	0.28,	0.0,  	0.28,  	0.28,],  # Row B
        [2.8,  	1.4,  	0.56,  	2.8,  	1.4,  	0.56,  	2.8,  	1.4,  	0.56,  	2.8,  	1.4,  	0.56,],  # Row C
        [0.028,	0.028,  2.8,    0.028,	0.028,  2.8,    0.028,	0.028,  2.8,   	0.028,	0.028,  2.8,],  # Row D
        [1.4,  	2.8,  	0.0,  	1.4,  	2.8,  	0.0, 	1.4,  	2.8,  	0.0,    1.4,  	2.8,  	0.0,],  # Row E
        [0.28,	0.0,    2.8,    0.28,	0.0,    2.8,  	0.28,	0.0,    2.8,  	0.28,	0.0,    2.8,],  # Row F
        [0.56,	2.8,    0.028,  0.56,	2.8,    0.028, 	0.56,	2.8,    0.028,  0.56,	2.8,    0.028,],  # Row Gh
        [2.8,	2.8,    0.0,    0.0,    2.8,    2.8,    0.0,    0.0,    2.8,    0.0,    0.0,    2.8],  # Row H
    ]),
	# Channel 2 (Color 2 or 6). Yellow-Green or White on v0.4c
	[[0 / (2**row) for col in range(NCOLS)] for row in range(NROWS)],    

    # Channel 3 – Red light (opposite checkerboard cells)
    np.array([
        [0.0,	2.8,	2.8,  	0.0,  	2.8,  	2.8,  	0.0,  	2.8,  	2.8,  	0.0,  	2.8,    2.8],  # Row A
        [2.8,   2.8,    2.8,    2.8,    2.8,    2.8,    2.8,    2.8,    2.8,   	2.8,    2.8,    2.8],  # Row B
        [2.8,   2.8,    2.8,    2.8,    2.8,    2.8,    2.8,    2.8,    2.8,   	2.8,    2.8,    2.8],  # Row C
        [2.8,  	2.8,  	0.0,  	2.8,  	2.8,  	0.0,  	2.8,  	2.8,  	0.0,  	2.8,  	2.8,    0.0],  # Row D
        [2.8,   2.8,    2.8,    2.8,   	2.8,    2.8,  	2.8,  	2.8,  	2.8,  	2.8,  	2.8,    2.8],  # Row E
        [2.8,  	2.8,  	2.8,  	2.8,  	2.8,  	2.8, 	2.8,    2.8,    2.8,    2.8,    2.8,    2.8],  # Row F
        [2.8,   0.0,    2.8,    2.8,    0.0,    2.8,    2.8,   	0.0,   	2.8,    2.8,    0.0,    2.8],  # Row G
        [0.0,   0.0,    2.8,   	2.8, 	0.0,  	0.0, 	2.8,  	2.8,  	0.0,  	2.8,  	2.8,   	0.0],  # Row H
    ])
])

#plate map for cell type
class Cells(Enum):
    media = 0
    JBL001 = 1
    YX001 = 2
    YX002 = 3

class Rows(Enum):
    A = 0
    B = 1
    C = 2
    D = 3
    E = 4
    F = 5
    G = 6
    H = 7

medium = {
    0: "LB",
    1: "WM-met+",
    2: "WM-met-",
}


cell_map = np.array([
        #1      2       3       4       5       6       7       8       9       10      11      12
        [2,	    2,	    2,  	2,  	2,  	2,  	3,  	3,  	3,  	3,  	3,      3],  # Row A
        [2,	    2,	    2,  	2,  	2,  	2,  	3,  	3,  	3,  	3,  	3,      3],  # Row B
        [2,	    2,	    2,  	2,  	2,  	2,  	3,  	3,  	3,  	3,  	3,      3],  # Row C
        [2,	    2,	    2,  	2,  	2,  	2,  	3,  	3,  	3,  	3,  	3,      3],  # Row D
        [2,	    2,	    2,  	2,  	2,  	2,  	3,  	3,  	3,  	3,  	3,      3],  # Row E
        [2,	    2,	    2,  	2,  	2,  	2,  	3,  	3,  	3,  	3,  	3,      3],  # Row F
        [2,	    2,	    2,  	2,  	2,  	2,  	3,  	3,  	3,  	3,  	3,      3],  # Row G
        [1,	    1,	    1,  	1,  	1,  	1,  	1,  	1,  	0,  	0,  	0,      0],  # Row H
])

media_map = np.array([
        #1      2       3       4       5       6       7       8       9       10      11      12
        [1,	    2,	    1,  	2,  	1,  	2,  	1,  	2,  	1,  	2,  	1,      2],  # Row A
        [1,	    2,	    1,  	2,  	1,  	2,  	1,  	2,  	1,  	2,  	1,      2],  # Row B
        [1,	    2,	    1,  	2,  	1,  	2,  	1,  	2,  	1,  	2,  	1,      2],  # Row C
        [1,	    2,	    1,  	2,  	1,  	2,  	1,  	2,  	1,  	2,  	1,      2],  # Row D
        [1,	    2,	    1,  	2,  	1,  	2,  	1,  	2,  	1,  	2,  	1,      2],  # Row E
        [1,	    2,	    1,  	2,  	1,  	2,  	1,  	2,  	1,  	2,  	1,      2],  # Row F
        [1,	    2,	    1,  	2,  	1,  	2,  	1,  	2,  	1,  	2,  	1,      2],  # Row G
        [1,	    1,	    1,  	1,  	2,  	2,  	2,  	2,  	1,  	2,  	1,      2],  # Row H
])

#turn light array into labels - need to code
plate_map_new = {key : np.nan for key in key_wells}

#placing cell names into arrays in the plate_map
for row in range(0,len(cell_map)):
    row_letter = Rows(row).name
    for index in range(0,len(cell_map[row])):
        plate_map[str(row_letter) + str(index+1)].append(Cells(cell_map[row][index]).name)

#doing the same for media
for row in range(0,len(media_map)):
    row_letter = Rows(row).name
    for index in range(0,len(media_map[row])):
        plate_map[str(row_letter) + str(index+1)].append(medium[media_map[row][index]])

#doing the same for the color intensities
green_array = DEFAULT_OPTICAL_POWER[1]
for row in range(0,len(green_array)):
    row_letter = Rows(row).name
    for index in range(0,len(green_array[row])):
        plate_map[str(row_letter) + str(index+1)].append(green_array[row][index])

red_array = DEFAULT_OPTICAL_POWER[3]
for row in range(0,len(red_array)):
    row_letter = Rows(row).name
    for index in range(0,len(red_array[row])):
        plate_map[str(row_letter) + str(index+1)].append(red_array[row][index])

#making a concatinated string name for easy accessing
for key, item in plate_map.items():
    plate_map[key].append("_".join(map(str, item)))





#Use the extracted_combined file from plate_reader_extraction script to get the right format
filepath = "26-09-10 YX002 test 1/26-09-10 YX002 test 1_diya_extracted_combined.csv"
DefaultConfig.save_file_path = filepath.split("/")[0]

sorted_data_df, summary_df = load_data(filepath)


# Color and style map
style_map = {
    "marker_style":{
        "o":[]
    },
    "marker_color":{},
    "line_style":{},
    "line_color":{},

}


#custom color and linestyles
#marker style and color
markerstyle_map = {"JBL001":"s",      
                    "YX001":"o",
                    "YX002":"s",
                    "media":"^",
}

markercolor_map = {"JBL001":"black",      
                    "YX001":"blue",
                    "YX002":"blue",
                    "media":"black",
}

#linestyle 
linestyle_map = {"WM-met+": "solid",
                 "WM-met-": "dashed",
}

plot_exclude = {
    "cells":[],
    "media":[],
    "green_intensity":["0.0","1.4","0.56","0.028"],
    "red_intensity":[],
}

#od600 - average
line_color_map = {2.8: (0.0,1.0,0.0),
             1.4: (0.5,0.5,0.0),
             0.56:(0.6,0.4,0.0),
             0.28:(0.0,1.0,0.0),
             0.028:(0.0,1.0,0.0),
             0.0:(1.0,0.0,0.0),
}


alpha_map = {"2.8": 1,
                 "2.52": 0.95,
                 "2.24": 0.9,
                 "1.96": 0.8,
                 "1.68": 0.7,
                 "1.4": 0.6,
                 "1.12":0.5,
                 "0.84":0.4,
                 "0.56":0.3,
                 "0.28":0.2,
                 "0":0.1,
}


plot_exclude = {
    "cells":[ "media", "YX001"],
    "media":[],
    #"green_intensity":[2.8,1.4,0.56,0.028],
    #"green_intensity":[2.8,1.4,0.56,0.028],
    "green_intensity":[],
    "red_intensity":[],
}
DefaultConfig.plot_exclude = plot_exclude
#DefaultConfig.line_color_map = line_color_map
plot_timecourse(sorted_data_df, summary_df, "GFP 395nm/OD600", "average", title_extra= "YX002", save_image = True)
#plot_timecourse(sorted_data_df, summary_df, "GFP 395nm/OD600", "all", title_extra= " ", save_image = False)


#plot_by_intensity(sorted_data_df, summary_df, "OD600", "average", 5, title_extra= "t12 YX compare all", save_image = True)
