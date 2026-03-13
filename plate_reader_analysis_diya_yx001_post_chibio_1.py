from copy import deepcopy
import os
from pathlib import Path
from datetime import datetime
from os import listdir, path
from enum import Enum
from dataclasses import dataclass
from typing import Dict, Any, Optional
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from collections import defaultdict

#classes
@dataclass
class PlotConfig:
    markerstyle_map: defaultdict
    markercolor_map: defaultdict
    linestyle_map: defaultdict
    line_color_map: defaultdict
    line_alpha_map: defaultdict
    plot_exclude: Dict[str, list]
    alpha_used: bool
    save_file_path: str

# Functions
def default_plot_config(save_file_path: str = ".") -> PlotConfig:
    return PlotConfig(
        markerstyle_map = defaultdict(lambda: "o", {"JBL001":"s", "JBL137":"o", "media":"^"}),
        markercolor_map = defaultdict(lambda: "blue", {"JBL001":"black", "JBL137":"blue", "media":"black"}),
        linestyle_map = defaultdict(lambda: "solid", {"WM-met+": "solid", "WM-met-": "dashed", "LB": "solid"}),
        line_color_map = {lambda: "black"},  # build from data later
        line_alpha_map = {lambda: 1.0},
        plot_exclude = {"cells":[], "media":[], "green_intensity":[], "red_intensity":[]},
        alpha_used = False,
        save_file_path = save_file_path
    )
DefaultConfig = default_plot_config()

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

def plot_timecourse(dataframe, y_data, plot_type, ylabel: str | None = None, title: str | None = None, title_extra: str = "",
                            xlabel = "Time (hrs)",
                            config: PlotConfig = DefaultConfig, save_image = False):
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
        index_name_cells = row["cells"]
        index_name_media = row["media"]
        index_name_green_intensity = row["green_intensity"]
        index_name_red_intensity = row["red_intensity"]

        if (index_name_cells in config.plot_exclude["cells"]
            or index_name_media in config.plot_exclude["media"]
            or index_name_green_intensity in config.plot_exclude["green_intensity"]
            or index_name_red_intensity in config.plot_exclude["red_intensity"]):
            continue

        if config.alpha_used is True and index_name_red_intensity != 0:
            alpha = config.line_alpha_map[index_name_red_intensity]
        else:
            alpha = 1

        if plot_type == "average":
            axs.errorbar(row["timepoints"], row[f"{y_data}_average"],
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
                axs.plot(row["timepoints"], repeat,
                        
                            color = config.line_color_map[index_name_green_intensity],
                            marker = config.markerstyle_map[index_name_cells],
                            markerfacecolor = config.markercolor_map[index_name_cells],
                            markeredgecolor = config.markercolor_map[index_name_cells], markersize = 3.0,
                            linestyle = config.linestyle_map[index_name_media], linewidth = 1.0,
                            alpha = alpha,
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
        bbox_to_anchor=(1.02, 0.45),
        borderaxespad=0.0,
    )
    axs.add_artist(leg2)

    leg3 = axs.legend(
        handles=intensity_handles,
        title="Green light intensity",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.25),
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



def plot_by_intensity(dataframe, y_data, plot_type, data_timearray_loc, ylabel: str | None = None, title: str | None = None,
                              title_extra: str = "",
                              xlabel = "Green light intensity %",
                                config: PlotConfig = DefaultConfig, save_image = False):
    """
    Plots y_data over intensity at a chosen time point
    
    :param dataframe: dataframe where the data is
    :param y_data: data to plot
    :param plot_type: select whether 'average' are plotted or 'all' data
    :param data_timearray_loc: chosen time point by index
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
        title = f"{ylabel} - by intensity - {plot_type} - {title_extra}"
    else:
        title = title + " - " + title_extra

    by_intensity_df = dataframe[["cells","media","green_intensity","red_intensity",f"{y_data}_average",f"{y_data}_std", f"{y_data}_raw_array"]].copy()
    by_intensity_df[f"{y_data}_average"] = by_intensity_df[f"{y_data}_average"].apply(lambda x: x[data_timearray_loc])
    by_intensity_df[f"{y_data}_std"] = by_intensity_df[f"{y_data}_std"].apply(lambda x: x[data_timearray_loc])
    by_intensity_df[f"{y_data}_raw_array"] = by_intensity_df[f"{y_data}_raw_array"].apply(lambda x: [array[data_timearray_loc] for array in x])
    by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
    by_intensity_df["red_intensity_percentage"] = by_intensity_df["red_intensity"].apply(lambda x: 100*x/2.8)

    colors = [(1, 0, 0), (1, 0, 0)]  # Red (1,0,0) to Green (0,1,0)
    cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

    fig, axs = plt.subplots()
    for celltype in set(by_intensity_df["cells"]):
        if celltype in config.plot_exclude["cells"]:
            continue
        data_select_by_cell = by_intensity_df.loc[by_intensity_df["cells"] == celltype]

        for media in set(data_select_by_cell["media"]):
            if media in config.plot_exclude["media"]:
                continue
            data_select_by_media = data_select_by_cell.loc[data_select_by_cell["media"] == media]

            if plot_type == "average":
                i_range = range(1)
            elif plot_type == "all":
                lengths = [len(r) for r in data_select_by_media[f"{y_data}_raw_array"]]
                max_len = int(max(lengths))
                i_range = range(max_len)

            for i in i_range:

                x = data_select_by_media["green_intensity_percentage"].values.copy()
                
                if plot_type == "average":
                    y = data_select_by_media[f"{y_data}_average"].values.copy()
                elif plot_type == "all":
                    y = np.array([ (row[i] if i < len(row) else np.nan) for row in data_select_by_media[f"{y_data}_raw_array"] ])

                yerr = data_select_by_media[f"{y_data}_std"].values.copy()

                # move final point (green=2.8 & red=0) to the right
                is_final = (
                    (data_select_by_media["green_intensity"] == 2.8) &
                    (data_select_by_media["red_intensity"] == 0)
                )

                x[is_final.values] = x.min() - 10  # move to the left
                
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
                lc.set_array(x)
                axs.add_collection(lc)

                # Add scatter points with gradient colors
                scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5,
                                    marker = config.markerstyle_map[celltype],
                                    facecolor = config.markercolor_map[celltype],
                                    edgecolor = config.markercolor_map[celltype],
                                    )

                # Add error bars
                if plot_type == "average":
                    axs.errorbar(x, y, yerr=yerr, fmt='none',
                                ecolor= config.markercolor_map[celltype], alpha=1.0, capsize=2.0,
                                )

    # labels
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
        bbox_to_anchor=(1.02, 0.45),
        borderaxespad=0.0,
    )
    axs.add_artist(leg2)

    x_axis = axs.set_xlabel(xlabel)
    x_axis.set_color("red")

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

#generating plate map
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

    # Channel 1 – Green light 
    np.array([
        [1.4,	1.4,	0.0,  	0.0,  	0.0,  	0.28,  	1.4,  	0.0,  	1.4,  	0.0,  	0.0,    1.4],  # Row A
        [0.28,	0.28,	0.28,  	0.28,  	1.4,  	1.4,  	0.0,  	0.28,  	0.28,  	1.4,  	0.28,    0.28],  # Row B
        [0.0,	0.0,	1.4,  	1.4,  	0.28,  	0.0,  	0.28,  	1.4,  	0,  	0.28,  	1.4,    0],  # Row C
        [1.4,	0.28,	0.28,  	1.4,  	1.4,  	0.28,  	1.4,  	0.0,  	1.4,  	0.28,  	0.0,    1.4],  # Row D
        [0.28,	1.4,	0.0,  	0.28,  	0.0,  	0.0,  	0.0,  	0.28,  	0.0,  	1.4,  	0.28,    0.0],  # Row E
        [0.0,	0.0,	1.4,  	0.0,  	1.4,  	1.4,  	0.0,  	1.4,  	0.28,  	0.0,  	1.4,    0.28],  # Row F
        [1.4,	0.28,	0.0,  	0.28,  	0.28,  	0.28,  	0.28,  	0.28,  	1.4,  	0.0,  	0.0,    0.0],  # Row G
        [0.0,  1.4,    0.28,   	0.0, 	0.0,  	1.4, 	1.4,  	0.0,  	0.28,  	0.0,  	0.0,   	0.0],  # Row H
    ]),
	# Channel 2 (Color 2 or 6). Yellow-Green or White on v0.4c
	[[0 / (2**row) for col in range(NCOLS)] for row in range(NROWS)],    

    # Channel 3 – Red light
    np.array([
        # Columns 1–12 (A–L), Rows A–H
        [2.8,	2.8,   2.8,    2.8,	2.8,   2.8,    2.8,	2.8,   2.8,   	2.8,	2.8,   2.8],  # Row A
        [2.8,	2.8,   2.8,    2.8,	2.8,   2.8,    2.8,	2.8,   2.8,   	2.8,	2.8,   2.8],  # Row B
        [2.8,	2.8,   2.8,    2.8,	2.8,   2.8,    2.8,	2.8,   2.8,   	2.8,	2.8,   2.8],  # Row C
        [2.8,	2.8,   2.8,    2.8,	2.8,   2.8,    2.8,	2.8,   2.8,   	2.8,	2.8,   2.8],  # Row D
        [2.8,	2.8,   2.8,    2.8,	2.8,   2.8,    2.8,	2.8,   2.8,   	2.8,	2.8,   2.8],  # Row E
        [2.8,	2.8,   2.8,    2.8,	2.8,   2.8,    2.8,	2.8,   2.8,   	2.8,	2.8,   2.8],  # Row F
        [2.8,	2.8,   2.8,    2.8,	2.8,   2.8,    2.8,	2.8,   2.8,   	0.0,	0.0,   0.0],  # Row G
        [2.8,	2.8,   2.8,    2.8,	2.8,   2.8,    2.8,	2.8,   2.8,   	0.0,	0.0,   0.0],  # Row H
    ])
])

#plate map for cell type
class Cells(Enum):
    media = 0
    JBL001 = 1
    YX001 = 2

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
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,      2],  # Row A
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,      2],  # Row B
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,      2],  # Row C
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,      2],  # Row D
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,      2],  # Row E
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,      2],  # Row F
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	1,  	1,      1],  # Row G
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	0,  	0,      0],  # Row H
])

media_map = np.array([
        #1      2       3       4       5       6       7       8       9       10      11      12
        [1,	    2,	    1,  	2,  	2,  	2,  	2,  	2,  	1,  	2,  	1,      2],  # Row A
        [1,	    2,	    1,  	2,  	2,  	2,  	2,  	2,  	1,  	2,  	1,      2],  # Row B
        [1,	    2,	    1,  	2,  	2,  	2,  	2,  	2,  	1,  	2,  	1,      2],  # Row C
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,      2],  # Row D
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,      2],  # Row E
        [2,	    2,	    2,  	2,  	1,  	2,  	1,  	2,  	2,  	2,  	2,      2],  # Row F
        [2,	    2,	    2,  	2,  	1,  	2,  	1,  	2,  	2,  	2,  	2,      1],  # Row G
        [2,     2,      2,   	2, 	    1,  	2, 	    1,  	2,  	2,  	2,  	2,   	2],  # Row H
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






#filepaths = ["25-11-19_diya_4/25-11-18_diya2_dose_curve_low_gain_extracted_OD600.csv","25-11-19_diya_4/25-11-18_diya2_dose_curve_low_gain_extracted_GFP 488nm.csv", "25-11-19_diya_4/25-11-18_diya2_dose_curve_low_gain_extracted_GFP 395nm.csv"]

filepaths = ["26-03-10_YX001_diya_then_diya/26-03-10_YX001_diya_then_diya_extracted_GFP 395nm.csv",
             "26-03-10_YX001_diya_then_diya/26-03-10_YX001_diya_then_diya_extracted_GFP 488nm.csv", 
             "26-03-10_YX001_diya_then_diya/26-03-10_YX001_diya_then_diya_extracted_OD600.csv"]


#initiating data df
indexes = sorted(set([value[-1] for value in plate_map.values()]))

sorted_data_df = pd.DataFrame(0.0, index=indexes, columns=["cells",
                                                           "media",
                                                           "green_intensity",
                                                           "red_intensity",
                                                           "timepoints",
                                                           "OD600_raw_array","OD600_average","OD600_std",
                                                           "GFP488_raw_array","GFP488_average","GFP488_std",
                                                           "GFP395_raw_array","GFP395_average","GFP395_std",
                                                           "GFP/OD600_raw_array","GFP/OD600_average","GFP/OD600_std"]
                                                            ).astype(object)
#initialising arrays
sorted_data_df["OD600_raw_array"] = [[] for _ in range(len(sorted_data_df))]
sorted_data_df["GFP488_raw_array"] = [[] for _ in range(len(sorted_data_df))]
sorted_data_df["GFP395_raw_array"] = [[] for _ in range(len(sorted_data_df))]
sorted_data_df["GFP/OD600_raw_array"] = [[] for _ in range(len(sorted_data_df))]

#reading extracted data file csv
file_list = {}
for filepath in filepaths:
    filename = filepath.split("_")[-1].replace(".csv","")
    data = pd.read_csv(filepath, header = 0, index_col= 0)
    
    #working out the time points
    timepoints = list(data.columns.values)
    try:
        initial_time = datetime.strptime(timepoints[0], "%Y-%m-%d %H:%M:%S")
    except ValueError:
        initial_time = datetime.strptime(timepoints[0], "%d/%m/%Y %H:%M:%S")

    timepoints_hrs = []
    for item in timepoints:
        try:
            timepoints_datetime = datetime.strptime(item, "%Y-%m-%d %H:%M:%S")
        except ValueError:
            timepoints_datetime = datetime.strptime(item, "%d/%m/%Y %H:%M:%S")
        time_delta = timepoints_datetime - initial_time
        timepoints_hrs.append(time_delta.total_seconds() / 3600)

    #renaming file column to the time
    data = data.rename(columns={old:new for (old,new) in zip(timepoints, timepoints_hrs)})
    
    file_list[filename] = data


#populate raw data from wells
for well_coord, value in plate_map.items():
    sorted_data_df.loc[value[-1],"cells"] = value[0]
    sorted_data_df.loc[value[-1],"media"] = value[1]
    sorted_data_df.loc[value[-1],"green_intensity"] = value[2]
    sorted_data_df.loc[value[-1],"red_intensity"] = value[3]
    sorted_data_df.loc[value[-1],"OD600_raw_array"].append(np.array(file_list["OD600"].loc[str(well_coord)]))
    sorted_data_df.loc[value[-1],"GFP488_raw_array"].append(np.array(file_list["GFP 488nm"].loc[str(well_coord)]))
    sorted_data_df.loc[value[-1],"GFP395_raw_array"].append(np.array(file_list["GFP 395nm"].loc[str(well_coord)]))
    sorted_data_df.loc[value[-1],"GFP/OD600_raw_array"].append(np.array(file_list["GFP 395nm"].loc[str(well_coord)])/np.array(file_list["OD600"].loc[str(well_coord)]))


#calculating time points
sorted_data_df["timepoints"] = [list(file_list["OD600"].columns.values) for _ in range(len(sorted_data_df))]

#calculating means and std
for index, row in sorted_data_df.iterrows():

    #mean, std data
    sorted_data_df.at[index,"OD600_average"] = np.mean(sorted_data_df.loc[index,"OD600_raw_array"], axis = 0)
    sorted_data_df.at[index,"OD600_std"] = np.std(sorted_data_df.loc[index,"OD600_raw_array"], axis = 0)
    sorted_data_df.at[index,"GFP488_average"] = np.mean(sorted_data_df.loc[index,"GFP488_raw_array"], axis = 0)
    sorted_data_df.at[index,"GFP488_std"] = np.std(sorted_data_df.loc[index,"GFP488_raw_array"], axis = 0)
    sorted_data_df.at[index,"GFP395_average"] = np.mean(sorted_data_df.loc[index,"GFP395_raw_array"], axis = 0)
    sorted_data_df.at[index,"GFP395_std"] = np.std(sorted_data_df.loc[index,"GFP395_raw_array"], axis = 0)

    sorted_data_df.at[index,"GFP/OD600_average"] = np.mean(sorted_data_df.loc[index,"GFP/OD600_raw_array"], axis = 0)
    sorted_data_df.at[index,"GFP/OD600_std"] = np.std(sorted_data_df.loc[index,"GFP/OD600_raw_array"], axis = 0)


#plotting

#default color and linestyles
#marker style and color
markerstyle_map = {"JBL001":"s",      
                    "YX001": "o",
                    "media":"^",
}


markercolor_map = {"JBL001":"black",      
                    "YX001":"blue",
                    "media":"black",
}

#linestyle 
linestyle_map = {"WM-met+": "solid",
                 "WM-met-": "dashed",
}

#default line_color is just percentage of green light -> (R,G,B)
line_color = lambda green_intensity: (1.0 - green_intensity/2.8, green_intensity/2.8, 0)
line_color_map = defaultdict(lambda:"black",{green_intensity : line_color(green_intensity) for green_intensity in set(sorted_data_df["green_intensity"])})

#default alpha is also percentage of green light, if used
#line_alpha = lambda green_intensity: green_intensity/2.8
#line_alpha_map = {green_intensity : line_alpha(green_intensity) for green_intensity in set(sorted_data_df["green_intensity"])}

#alpha using red intensity instead
line_alpha = lambda green_intensity: green_intensity/2.8
line_alpha_map = {green_intensity : line_alpha(green_intensity) for green_intensity in set(sorted_data_df["green_intensity"])}


DefaultConfig.line_color_map = line_color_map
DefaultConfig.line_alpha_map = line_alpha_map

plot_exclude = {
    "cells":[],
    "media":[],
    "green_intensity":[],
    "red_intensity":[],
}

#od600 - average
line_color_map_override = {"2.8": (0.0,1.0,0.0),
             "1.4": (0.5,0.5,0.0),
             "0.56":(0.6,0.4,0.0),
             "0.28":(0.8,0.2,0.0),
             "0.028":(0.9,0.1,0.0),
             "0.0":(1.0,0.0,0.0),
}


alpha_map_override = {"2.8": 1,
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

save_file_path = "26-03-10_YX001_diya_then_diya"

#legend handles
cell_handles = [
    Line2D(
        [], [], 
        marker = markerstyle_map[cell],
        linestyle = "None",
        markerfacecolor = markercolor_map[cell],
        markeredgecolor = markercolor_map[cell],
        markersize = 6,
        label = cell
    )
    for cell in markerstyle_map
]

sorted_intensities = sorted(
    line_color_map.keys(),
    key=float,
    reverse=True
)

intensity_handles = [
    Line2D(
        [], [],
        color=line_color_map[intensity],
        linewidth=2,
        label=intensity
    )
    for intensity in sorted_intensities
]

media_handles = [
    Line2D(
        [], [],
        color = "black",
        linestyle = linestyle_map[media],
        linewidth = 2,
        label = media
    )
    for media in linestyle_map
]

alpha_handles = [
    Line2D(
        [], [],
        color = line_color_map[intensity],
        alpha = line_alpha_map[intensity],
        linewidth=2,
        label=intensity
    )
    for intensity in sorted_intensities
]







plot_exclude_override = {
    "cells":[
            #"YX001",
            "JBL001",
            ],
    "media":[
            #"WM-met-",
            "WM-met+",
        ],
    "green_intensity":[],
    "red_intensity":[],
}

DefaultConfig.save_file_path = "26-03-10_YX001_diya_then_diya"
OverrideConfig = deepcopy(DefaultConfig)
#OverrideConfig.alpha_used = True
OverrideConfig.markerstyle_map = markerstyle_map
OverrideConfig.markercolor_map = markercolor_map
OverrideConfig.plot_exclude = plot_exclude_override

for i in ["OD600", "GFP395", "GFP488", "GFP/OD600"]:
    #_title_extra = "YX001 post chibio 2"
    _title_extra = "met -"
    plot_timecourse(sorted_data_df, i, "average", config = OverrideConfig, title_extra= _title_extra, save_image = True)
    plot_timecourse(sorted_data_df, i, "all", config = OverrideConfig, title_extra= _title_extra, save_image = True)

    #_title_extra = "YX001 post chibio 2, t12"
    #plot_by_intensity(sorted_data_df,i, "average", -2, config = OverrideConfig, title_extra= _title_extra, save_image = True)
    #plot_by_intensity(sorted_data_df,i, "all", -2, config = OverrideConfig, title_extra= _title_extra, save_image = True)
    pass




def plot_by_intensity_all_separate(dataframe, y_data, data_timearray_loc, plot_select,
                                   #plot_control,
                                   ylabel: str | None = None, title: str | None = None, title_extra: str = "", xlabel = "Red light intensity %",
                        plot_exclude = plot_exclude, markerstyle_map = markerstyle_map,
                        linestyle_map = linestyle_map, line_color_map = line_color_map, line_alpha_map = line_alpha_map,
                        alpha_used = False, save_image = False, save_filepath = None):
    
    if ylabel is None:
        ylabel = y_data

    if title is None:
        title = f"{ylabel} - by intensity - all - {title_extra}"
    else:
        title = title + " - " + title_extra

    if save_filepath is None:
        save_filepath = save_file_path

    by_intensity_df = dataframe[["cells","media","green_intensity","red_intensity",f"{y_data}_raw_array"]].copy()
    by_intensity_df[f"{y_data}_raw_array"] = by_intensity_df[f"{y_data}_raw_array"].apply(lambda x: [array[data_timearray_loc] for array in x])
    by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
    by_intensity_df["red_intensity_percentage"] = by_intensity_df["red_intensity"].apply(lambda x: 100*x/2.8)

    colors = [(1, 0, 0), (1, 0, 0)]  # Red (1,0,0) to Green (0,1,0)
    cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)


    data_select = by_intensity_df[by_intensity_df["cells"].isin(plot_select["cells"]) & by_intensity_df["media"].isin(plot_select["media"])]
    n_positions = len(data_select[f"{y_data}_raw_array"].iloc[0]) if len(data_select) > 0 else 0
    
    #data_control = by_intensity_df[by_intensity_df["cells"].isin(plot_control["cells"]) and by_intensity_df["media"].isin(plot_control["media"])]

    n_cols = 2  # 2 columns
    n_rows = (n_positions + n_cols - 1) // n_cols  # Calculate rows needed
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))
    axs = axs.flatten()  # Flatten to 1D array for easier indexing

    for i in range(n_positions):
        ax = axs[i]
        x = data_select["red_intensity_percentage"].values.copy()
        y = np.array([ (row[i] if i < len(row) else np.nan) for row in data_select[f"{y_data}_raw_array"] ])

        # move final point (green=2.8 & red=0) to the right
        is_final = (
            (data_select["green_intensity"] == 2.8) &
            (data_select["red_intensity"] == 0)
        )

        x[is_final.values] = x.min() - 10  # move to the right
        
        order = np.argsort(x)
        x = x[order]
        y = y[order]

        # Create line segments for gradient
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Normalize x values for colormap
        media = plot_select["media"][0]
        celltype = plot_select["cells"][0]

        norm = plt.Normalize(x.min(), x.max())
        lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
        lc.set_linestyle(linestyle_map[media])
        lc.set_array(x)
        ax.add_collection(lc)

        # Add scatter points with gradient colors
        scatter = ax.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5,
                            marker = markerstyle_map[celltype],
                            facecolor = markercolor_map[celltype],
                            edgecolor = markercolor_map[celltype],
                            )
        
        # Auto-scale to fit LineCollection
        ax.autoscale()
        ax.set_title(f"Repeat {i+1}", fontsize = 8)
        ax.grid(True, alpha=0.3)

    # labels
    leg1 = ax.legend(
        handles=cell_handles,
        title="Cell type",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.10),
        borderaxespad=0.0,
    )
    ax.add_artist(leg1)

    leg2 = ax.legend(
        handles=media_handles,
        title="Media",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.65),
        borderaxespad=0.0,
    )
    ax.add_artist(leg2)

    fig.supxlabel(xlabel)
    fig.supylabel(ylabel)
    fig.suptitle(f"{ylabel} by intensity - all repeats separated", fontsize=14, y=1.00)

    # Hide any unused subplots
    for j in range(n_positions, len(axs)):
        axs[j].axis('off')

    fig.tight_layout()

    if save_image == True:
        try:
            Path(os.path.join(save_filepath, "figures")).mkdir(parents = True, exist_ok = True)
            save_title = title.replace("/","_div_")
            plt.savefig(os.path.join(save_filepath, "figures", f"{save_title}.png"))
            plt.close(fig)
        except:
            print("Could not save image, filepath not valid")
    else:
        plt.show()
        plt.close(fig)


plot_select_override = {"cells":["JBL137"],"media":["WM-met-"]}
#plot_by_intensity_all_separate(sorted_data_df, "GFP395", -2, plot_select_override, title_extra="WM-met- t12", save_image=True, save_filepath="26-01-15_new_jbl137_diya_wm_red")
#plot_by_intensity_all_separate(sorted_data_df, "OD600", -2, plot_select_override, title_extra="WM-met- t12", save_image=True, save_filepath="26-01-15_new_jbl137_diya_wm_red")
#plot_by_intensity_all_separate(sorted_data_df, "GFP/OD600", -2, plot_select_override, title_extra="WM-met- t12", save_image=True, save_filepath="26-01-15_new_jbl137_diya_wm_red")
plot_select_override = {"cells":["JBL137"],"media":["WM-met+"]}
#plot_by_intensity_all_separate(sorted_data_df, "GFP395", -2, plot_select_override, title_extra="WM-met+ t12", save_image=True, save_filepath="26-01-15_new_jbl137_diya_wm_red")
#plot_by_intensity_all_separate(sorted_data_df, "OD600", -2, plot_select_override, title_extra="WM-met+ t12", save_image=True, save_filepath="26-01-15_new_jbl137_diya_wm_red")
#plot_by_intensity_all_separate(sorted_data_df, "GFP/OD600", -2, plot_select_override, title_extra="WM-met+ t12", save_image=True, save_filepath="26-01-15_new_jbl137_diya_wm_red")




####
#GFP by intensity - average over time
def plot_by_intensity_average_over_time(dataframe, y_data, plot_select, ylabel: str | None = None, title: str | None = None, title_extra: str = "", xlabel = "Red light intensity %",
                        plot_exclude = plot_exclude, markerstyle_map = markerstyle_map,
                        linestyle_map = linestyle_map, line_color_map = line_color_map, line_alpha_map = line_alpha_map,
                        alpha_used = False, save_image = False, save_filepath = None):
    
    if ylabel is None:
        ylabel = y_data

    if title is None:
        title = f"{ylabel} - by intensity - averages over time - {title_extra}"
    else:
        title = title + " - " + title_extra

    if save_filepath is None:
        save_filepath = save_file_path

    by_intensity_df = dataframe[["cells","media","timepoints","green_intensity","red_intensity",f"{y_data}_average",f"{y_data}_std"]].copy()
    by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
    by_intensity_df["red_intensity_percentage"] = by_intensity_df["red_intensity"].apply(lambda x: 100*x/2.8)

    times = by_intensity_df["timepoints"][0]
    times_rounded = [round(x) for x in times]

    err_norm = mcolors.Normalize(vmin=min(times), vmax=max(times))
    err_cmap = plt.cm.tab20c
    colors = [(1, 0, 0), (1, 0, 0)]  # Red (1,0,0) to Green (0,1,0)
    cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

    errorbar_handles = []
    for t in times_rounded:
        color = err_cmap(err_norm(t))
        handle = Line2D(
            [0], [0],
            color=color,
            linewidth=2,
            linestyle='-',
            label=f"t{t}"
            )
        errorbar_handles.append(handle)

    data_select = by_intensity_df[by_intensity_df["cells"].isin(plot_select["cells"]) & by_intensity_df["media"].isin(plot_select["media"])]
    fig, axs = plt.subplots()

    for celltype in set(data_select["cells"]):
        data_select_by_cell = data_select.loc[data_select["cells"] == celltype]

        for media in set(data_select_by_cell["media"]):
            data_select_by_media = data_select_by_cell.loc[data_select_by_cell["media"] == media]

            for i in range(len(times_rounded)):
                data_select_by_time_avg = data_select_by_media[f"{y_data}_average"].copy().apply(lambda x: x[i])
                data_select_by_time_std = data_select_by_media[f"{y_data}_std"].copy().apply(lambda x: x[i])
                # selecting data based on time

                # plot with gradient color
                x = data_select_by_media["red_intensity_percentage"].values.copy()
                y = data_select_by_time_avg
                yerr = data_select_by_time_std

                # move final point (green=2.8 & red=0) to the right
                is_final = (
                    (data_select_by_media["green_intensity"] == 2.8) &
                    (data_select_by_media["red_intensity"] == 0)
                )

                x[is_final.values] = x.min() - 10  # move to the right
                
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
                lc.set_linestyle(linestyle_map[media])
                lc.set_array(x)
                axs.add_collection(lc)

                # Add scatter points with gradient colors
                scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5,
                                    marker = markerstyle_map[celltype],
                                    facecolor = markercolor_map[celltype],
                                    edgecolor = markercolor_map[celltype],
                                    )

                # Add error bars
                err_color = err_cmap(err_norm(times_rounded[i]))
                axs.errorbar(x, y, yerr=yerr, fmt='none',
                            ecolor= err_color, alpha=1.0, capsize=2.0,
                            label = "t"+ str(times_rounded[i]),
                            )

    # labels
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
        handles=errorbar_handles,
        title="Time (error bars)",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.30),
        borderaxespad=0.0,
    )
    axs.add_artist(leg3)

    x_axis = axs.set_xlabel(xlabel)
    x_axis.set_color("red")

    axs.set_ylabel(ylabel)
    axs.set_title(title)

    fig.tight_layout(rect=[0, 0, 0.75, 1])
    if save_image == True:
        try:
            Path(os.path.join(save_filepath, "figures")).mkdir(parents = True, exist_ok = True)
            save_title = title.replace("/","_div_")
            plt.savefig(os.path.join(save_filepath, "figures", f"{save_title}.png"))
            plt.close(fig)
        except:
            print("Could not save image, filepath not valid")
    else:
        plt.show()
        plt.close(fig)


plot_select_override = {"cells":["JBL137"],"media":["WM-met-","WM-met+"]}
#plot_by_intensity_average_over_time(sorted_data_df, "OD600", plot_select_override, title_extra="all media", save_image=True, save_filepath="26-01-15_new_jbl137_diya_wm_red")






def plot_by_intensity_foldchange(dataframe, y_data, data_timearray_loc, ylabel: str | None = None, title: str | None = None, title_extra: str = "", xlabel = "Red light intensity %",
                        plot_exclude = plot_exclude, markerstyle_map = markerstyle_map,
                        linestyle_map = linestyle_map, line_color_map = line_color_map, line_alpha_map = line_alpha_map,
                        alpha_used = False, save_image = False, save_filepath = None):
#calculating fold changes
    if ylabel is None:
        ylabel = y_data

    if title is None:
        title = f"{ylabel} - by intensity - fold change - {title_extra}"
    else:
        title = title + " - " + title_extra

    if save_filepath is None:
        save_filepath = save_file_path



    data_select = dataframe[["cells","media", "green_intensity", "red_intensity",f"{y_data}_average",f"{y_data}_std"]]
    data_select["avg_select"] = data_select[f"{y_data}_average"].copy().apply(lambda x: x[data_timearray_loc])
    data_select["std_select"] = data_select[f"{y_data}_std"].copy().apply(lambda x: x[data_timearray_loc])

    media_row = data_select.loc[data_select["cells"] == "media"]
    media_value = media_row[f"{y_data}_average"].copy().apply(lambda x: x[data_timearray_loc])

    data_select["avg_normalised"] = data_select["avg_select"] - media_value.values
    data_select = data_select[data_select["cells"] != "media"]

    met_plus_df = data_select[data_select["media"] == "WM-met+"].copy()
    met_minus_df = data_select[data_select["media"] == "WM-met-"].copy()

    met_plus_df = met_plus_df[["cells", "green_intensity", "red_intensity", "avg_normalised", "std_select"]].rename(
        columns={"avg_normalised": "avg_met_plus", "std_select": "std_met_plus"}
    )
    met_minus_df = met_minus_df[["cells", "green_intensity", "red_intensity", "avg_normalised", "std_select"]].rename(
        columns={"avg_normalised": "avg_met_minus", "std_select": "std_met_minus"}
    )

    paired = pd.merge(
        met_minus_df,
        met_plus_df,
        on=["cells", "green_intensity", "red_intensity"],
        how="inner"  # use 'outer' to keep unpaired rows if you want to inspect them
    )

    paired["fold_change"] = paired["avg_met_plus"] / paired["avg_met_minus"]
    paired["fc_std"] = paired["fold_change"] * np.sqrt(
        (paired["std_met_minus"] / paired["avg_met_minus"].replace(0, np.nan))**2 +
        (paired["std_met_plus"]  / paired["avg_met_plus"])**2
    )

    paired["red_intensity_percentage"] = paired["red_intensity"].apply(lambda x: 100*x/2.8)

    cells = paired["cells"].unique()
    palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, ax = plt.subplots(figsize=(7,5))

    for ci, cell in enumerate(cells):
        dfc = paired[paired["cells"] == cell].copy()
        if dfc.empty:
            continue

        # base x
        x = dfc["red_intensity_percentage"].to_numpy(dtype=float)

        # shift where red intensity == 0 by +10
        final_mask = dfc["green_intensity"].to_numpy(dtype=float) == 0
        x_shift = x.copy()
        x_shift[final_mask] = x_shift[final_mask] - 10.0

        # y and yerr
        y = dfc["fold_change"].to_numpy(dtype=float)
        yerr = dfc["fc_std"].to_numpy(dtype=float)

        # sort by x_shift so lines look correct left->right
        order = np.argsort(x_shift)
        x_s = x_shift[order]
        y_s = y[order]
        yerr_s = yerr[order]

        # line (use one palette color per cell)
        line_color = palette[ci % len(palette)]
        ax.plot(x_s, y_s, linestyle='-', marker=None, color=line_color, label=str(cell))

        # scatter markers: use markerstyle_map and markercolor_map
        marker = markerstyle_map.get(cell, "o")
        mcolor = markercolor_map.get(cell, line_color)
        ax.scatter(x_s, y_s, marker=marker, color=mcolor, edgecolors='k', zorder=5, s=50)

        # errorbars
        ax.errorbar(x_s, y_s, yerr=yerr_s, fmt='none', ecolor=mcolor, alpha=0.9, capsize=3)


    # Add legends for cell types (markers)
    cell_handles = [
        Line2D([0],[0], marker=markerstyle_map.get(c, "o"), color="w",
            markerfacecolor=markercolor_map.get(c, palette[i % len(palette)]),
            markeredgecolor='k', markersize=7, label=str(c))
        for i,c in enumerate(cells)
    ]
    leg1 = ax.legend(handles=cell_handles, title="Cell type", loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.add_artist(leg1)

    # axes labels / title / grid
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Fold change (met+ / met-)")
    ax.set_title(title)
    ax.grid(alpha=0.25)

    # tighten layout and show
    fig.tight_layout(rect=[0,0,0.78,1])  # leave room for legends on right
    if save_image == True:
        try:
            Path(os.path.join(save_filepath, "figures")).mkdir(parents = True, exist_ok = True)
            save_title = title.replace("/","_div_")
            plt.savefig(os.path.join(save_filepath, "figures", f"{save_title}.png"))
            plt.close(fig)
        except:
            print("Could not save image, filepath not valid")
    else:
        plt.show()
        plt.close(fig)


plot_by_intensity_foldchange(sorted_data_df, "OD600", -2, save_image=True, save_filepath="26-01-15_new_jbl137_diya_wm_red")