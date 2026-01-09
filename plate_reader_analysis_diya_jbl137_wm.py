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
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from enum import Enum

#functions and defaults
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
        [2.8,	2.8,    0.0,    0.0,    0.0,    2.8,    2.8,    0.0,    0.0,    0.0,    0.0,    0.0],  # Row H
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
        [0.0,   0.0,    2.8,   	2.8, 	0.0,  	0.0, 	0.0,  	2.8,  	2.8,  	0.0,  	2.8,   	2.8],  # Row H
    ])
])

#plate map for cell type
class Cells(Enum):
    media = 0
    JBL001 = 1
    JBL137 = 2

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
        [2,	    2,	    2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,  	2,      2],  # Row G
        [1,     1,      1,   	1, 	    1,  	1, 	    1,  	1,  	1,  	1,  	0,   	0],  # Row H
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
        [1,     1,      1,   	1, 	    1,  	1, 	    2,  	2,  	2,  	2,  	2,   	2],  # Row H
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

filepaths = ["26-01-07_new_jbl137_diya_wm/26-01-07_wm_extracted_GFP 395nm.csv",
             "26-01-07_new_jbl137_diya_wm/26-01-07_wm_extracted_GFP 488nm.csv", 
             "26-01-07_new_jbl137_diya_wm/26-01-07_wm_extracted_OD600.csv"]


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
                    "JBL137":"o",
                    "media":"^",
}

markercolor_map = {"JBL001":"black",      
                    "JBL137":"blue",
                    "media":"black",
}

#linestyle 
linestyle_map = {"WM-met+": "solid",
                 "WM-met-": "dashed",
}

#default line_color is just percentage of green light -> (R,G,B)
line_color = lambda green_intensity: (1.0 - green_intensity/2.8, green_intensity/2.8, 0)
line_color_map = {green_intensity : line_color(green_intensity) for green_intensity in set(sorted_data_df["green_intensity"])}

#default alpha is also percentage of green light, if used
line_alpha = lambda green_intensity: green_intensity/2.8
line_alpha_map = {green_intensity : line_alpha(green_intensity) for green_intensity in set(sorted_data_df["green_intensity"])}


plot_exclude = {
    "cells":["JBL001","media"],
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

save_file_path = "26-01-07_new_jbl137_diya_wm"

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



def plot_average_timecourse(dataframe, y_data, ylabel: str | None = None, title: str | None = None, xlabel = "Time (hrs)",
                            plot_exclude = plot_exclude, markerstyle_map = markerstyle_map, markercolor_map = markercolor_map,
                            linestyle_map = linestyle_map, line_color_map = line_color_map, line_alpha_map = line_alpha_map,
                            alpha_used = False, save_image = False, save_filepath = None):

    if ylabel is None:
        ylabel = y_data

    if title is None:
        title = f"{ylabel} - average"

    if alpha_used == True:
        alpha = line_alpha_map[index_name_green_intensity]
    else:
        alpha = 1

    if save_filepath is None:
        save_filepath = save_file_path

    fig, axs = plt.subplots()
    for row_index, row in dataframe.iterrows():
        index_name_cells = row["cells"]
        index_name_media = row["media"]
        index_name_green_intensity = row["green_intensity"]
        index_name_red_intensity = row["red_intensity"]
        if index_name_cells in plot_exclude["cells"] or index_name_media in plot_exclude["media"] or index_name_green_intensity in plot_exclude["green_intensity"]:
            continue
        else:
            axs.errorbar(row["timepoints"], row[f"{y_data}_average"],
                        yerr = row[f"{y_data}_std"], capsize = 2.0,
                        #label = row_index,
                        color = line_color_map[index_name_green_intensity],
                        marker = markerstyle_map[index_name_cells],
                        markerfacecolor = markercolor_map[index_name_cells],
                        markeredgecolor = markercolor_map[index_name_cells], markersize = 3.0,
                        linestyle = linestyle_map[index_name_media], linewidth = 1.0,
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

    axs.set_xlabel("Time (hrs)")
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


def plot_all_timecourse(dataframe, y_data, ylabel: str | None = None, title: str | None = None, xlabel = "Time (hrs)",
                            plot_exclude = plot_exclude, markerstyle_map = markerstyle_map,
                            linestyle_map = linestyle_map, line_color_map = line_color_map, line_alpha_map = line_alpha_map,
                            alpha_used = False, save_image = False, save_filepath = None):

    if ylabel is None:
        ylabel = y_data

    if title is None:
        title = f"{ylabel} - all"

    if alpha_used == True:
        alpha = line_alpha_map[index_name_green_intensity]
    else:
        alpha = 1

    if save_filepath is None:
        save_filepath = save_file_path

    fig, axs = plt.subplots()
    for row_index, row in dataframe.iterrows():
        index_name_cells = row["cells"]
        index_name_media = row["media"]
        index_name_green_intensity = row["green_intensity"]
        index_name_red_intensity = row["red_intensity"]
        if index_name_cells in plot_exclude["cells"] or index_name_media in plot_exclude["media"] or index_name_green_intensity in plot_exclude["green_intensity"]:
            continue
        else:
            for repeat in row[f"{y_data}_raw_array"]:
                axs.plot(row["timepoints"], repeat,
                            color = line_color_map[index_name_green_intensity],
                            marker = markerstyle_map[index_name_cells],
                            markerfacecolor = markercolor_map[index_name_cells],
                            markeredgecolor = markercolor_map[index_name_cells], markersize = 3.0,
                            linestyle = linestyle_map[index_name_media], linewidth = 1.0,
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

    axs.set_xlabel("Time (hrs)")
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


plot_average_timecourse(sorted_data_df, "OD600", save_image = True)
plot_average_timecourse(sorted_data_df, "GFP395", save_image = True)
plot_average_timecourse(sorted_data_df, "GFP/OD600", save_image = True)
plot_average_timecourse(sorted_data_df, "GFP488", save_image = True)
plot_all_timecourse(sorted_data_df, "OD600", save_image = True)
plot_all_timecourse(sorted_data_df, "GFP395", save_image = True)
plot_all_timecourse(sorted_data_df, "GFP/OD600", save_image = True)
plot_all_timecourse(sorted_data_df, "GFP488", save_image = True)





#GFP by intensity - average

by_intensity_df = sorted_data_df[["green_intensity","GFP395_average","GFP395_std"]].copy()
by_intensity_df["GFP395_average"] = by_intensity_df["GFP395_average"].apply(lambda x: x[-1])
by_intensity_df["GFP395_std"] = by_intensity_df["GFP395_std"].apply(lambda x: x[-1])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jbl_new_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137new')]
jbl_hybrid_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137hybrid')]
jbl_kirill_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137kirill')]
jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
media_data = by_intensity_df[by_intensity_df.index.str.contains('media')]

colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

fig, axs = plt.subplots()


""" axs.errorbar(jbl_new_data["green_intensity_percentage"], jbl_new_data["GFP395_average"],
            yerr = jbl_new_data["GFP395_std"], capsize = 2.0,
            label = "JBL001",
            color = "black",
            marker = "^", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            ) """

# JCCO plot with gradient color
x = jbl_kirill_data["green_intensity_percentage"].values
y = jbl_kirill_data["GFP395_average"].values
yerr = jbl_kirill_data["GFP395_std"].values

# Create line segments for gradient
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Normalize x values for colormap
norm = plt.Normalize(x.min(), x.max())
lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
lc.set_array(x)
axs.add_collection(lc)

# Add scatter points with gradient colors
scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5, marker='o')

# Add error bars for JCCO
axs.errorbar(x, y, yerr=yerr, fmt='none', ecolor='gray', alpha=1.0, capsize=2.0)

""" axs.errorbar(media_data["green_intensity_percentage"], media_data["GFP395_average"],
            yerr = media_data["GFP395_std"], capsize = 2.0,
            label = "media",
            color = "black",
            marker = "^", markersize = 3.0,
            linestyle = "dotted", linewidth = 1.0,
            ) """

x_axis = axs.set_xlabel("Green light intensity %")
x_axis.set_color("green")

axs.set_ylabel("GFP 395nm")
axs.set_title(f"GFP 395nm by intensity - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()



#GFP/OD600 by intensity - average

by_intensity_df = sorted_data_df[["green_intensity","GFP/OD600_average","GFP/OD600_std"]].copy()
by_intensity_df["GFP/OD600_average"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[-1])
by_intensity_df["GFP/OD600_std"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[-1])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jbl_new_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137new')]
jbl_hybrid_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137hybrid')]
jbl_kirill_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137kirill')]
jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
media_data = by_intensity_df[by_intensity_df.index.str.contains('media')]

colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

fig, axs = plt.subplots()


""" axs.errorbar(jbl_data["green_intensity_percentage"], jbl_data["GFP/OD600_average"],
            yerr = jbl_data["GFP/OD600_std"], capsize = 2.0,
            label = "JBL001",
            color = "black",
            marker = "^", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            ) """

""" axs.errorbar(jcco_data["green_intensity_percentage"], jcco_data["GFP/OD600_average"],
            yerr = jcco_data["GFP/OD600_std"], capsize = 2.0,
            label = "JCCO",
            #color = "blue",
            marker = "o", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            cmap = cmap,
            ) """

# JCCO plot with gradient color
x = jbl_kirill_data["green_intensity_percentage"].values
y = jbl_kirill_data["GFP/OD600_average"].values
yerr = jbl_kirill_data["GFP/OD600_std"].values

# Create line segments for gradient
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Normalize x values for colormap
norm = plt.Normalize(x.min(), x.max())
lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
lc.set_array(x)
axs.add_collection(lc)

# Add scatter points with gradient colors
scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5, marker='o')

# Add error bars for JCCO
axs.errorbar(x, y, yerr=yerr, fmt='none', ecolor='gray', alpha=1.0, capsize=2.0)

""" axs.errorbar(media_data["green_intensity_percentage"], media_data["GFP/OD600_average"],
            yerr = media_data["GFP/OD600_std"], capsize = 2.0,
            label = "media",
            color = "black",
            marker = "^", markersize = 3.0,
            linestyle = "dotted", linewidth = 1.0,
            )
 """
x_axis = axs.set_xlabel("Green light intensity %")
x_axis.set_color("green")

axs.set_ylabel("GFP 395nm/OD600")
axs.set_title(f"GFP 395nm/OD600 by intensity - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()


#OD600 by intensity - average

by_intensity_df = sorted_data_df[["green_intensity","OD600_average","OD600_std"]].copy()
by_intensity_df["OD600_average"] = by_intensity_df["OD600_average"].apply(lambda x: x[-1])
by_intensity_df["OD600_std"] = by_intensity_df["OD600_std"].apply(lambda x: x[-1])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jbl_new_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137new')]
jbl_hybrid_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137hybrid')]
jbl_kirill_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137kirill')]
jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
media_data = by_intensity_df[by_intensity_df.index.str.contains('media')]

colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

fig, axs = plt.subplots()

# JCCO plot with gradient color
x = jbl_kirill_data["green_intensity_percentage"].values
y = jbl_kirill_data["OD600_average"].values
yerr = jbl_kirill_data["OD600_std"].values

# Create line segments for gradient
points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Normalize x values for colormap
norm = plt.Normalize(x.min(), x.max())
lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
lc.set_array(x)
axs.add_collection(lc)

# Add scatter points with gradient colors
scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5, marker='o')

# Add error bars for JCCO
axs.errorbar(x, y, yerr=yerr, fmt='none', ecolor='gray', alpha=1.0, capsize=2.0)


x_axis = axs.set_xlabel("Green light intensity %")
x_axis.set_color("green")

axs.set_ylabel("OD600")
axs.set_title(f"OD600 by intensity - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#GFP by intensity - all
by_intensity_df = sorted_data_df[["green_intensity","GFP395_raw_array"]].copy()
by_intensity_df["GFP395_raw_array"] = by_intensity_df["GFP395_raw_array"].apply(lambda x: [array[-1] for array in x])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jbl_new_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137new')]
jbl_hybrid_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137hybrid')]
jbl_kirill_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137kirill')]
jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
media_data = by_intensity_df[by_intensity_df.index.str.contains('media')]

colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

fig, axs = plt.subplots()

for i in range(len(jbl_hybrid_data["GFP395_raw_array"].iloc[0])):

    # JCCO plot with gradient color
    x = jbl_hybrid_data["green_intensity_percentage"].values
    y = np.array([row[i] for row in jbl_hybrid_data['GFP395_raw_array']])

    # Create line segments for gradient
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Normalize x values for colormap
    norm = plt.Normalize(x.min(), x.max())
    lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
    lc.set_array(x)
    axs.add_collection(lc)

    # Add scatter points with gradient colors
    label_jcco = f'jbl_new_data' if i == 0 else None   
    scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5, marker='o')

""" for i in range(len(jbl_data["GFP395_raw_array"].iloc[0])):
    axs.plot(jbl_data["green_intensity_percentage"], [row[i] for row in jbl_data['GFP395_raw_array']],
                label = "JBL001" if i == 0 else None,
                color = "black",
                marker = "^", markersize = 3.0,
                linestyle = "solid", linewidth = 1.0,
                ) """

x_axis = axs.set_xlabel("Green light intensity %")
x_axis.set_color("green")

axs.set_ylabel("GFP 395nm")
axs.set_title(f"GFP 395nm by intensity - all")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()


#GFP by intensity - all positions in separate subplots

by_intensity_df = sorted_data_df[["green_intensity","GFP395_raw_array"]].copy()
by_intensity_df["GFP395_raw_array"] = by_intensity_df["GFP395_raw_array"].apply(lambda x: [array[-1] for array in x])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jbl_new_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137new')]
jbl_hybrid_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137hybrid')]
jbl_kirill_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137kirill')]
jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
media_data = by_intensity_df[by_intensity_df.index.str.contains('media')]

colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

# Get max number of positions across all datasets
n_jcco = len(jbl_kirill_data["GFP395_raw_array"].iloc[0]) if len(jbl_kirill_data) > 0 else 0
#n_jbl = len(jbl_data["GFP395_raw_array"].iloc[0]) if len(jbl_data) > 0 else 0
#n_media = len(media_data["GFP395_raw_array"].iloc[0]) if len(media_data) > 0 else 0
n_positions = n_jcco

# Create subplots - adjust rows and cols as needed
n_cols = 2  # 2 columns
n_rows = (n_positions + n_cols - 1) // n_cols  # Calculate rows needed
fig, axs = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))
axs = axs.flatten()  # Flatten to 1D array for easier indexing

for i in range(n_positions):
    ax = axs[i]
    
    # JCCO plot with gradient color
    if i < n_jcco:
        x = jbl_kirill_data["green_intensity_percentage"].values
        y = np.array([row[i] for row in jbl_kirill_data['GFP395_raw_array']])

        # Create line segments for gradient
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Normalize x values for colormap
        norm = plt.Normalize(x.min(), x.max())
        lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
        lc.set_array(x)
        ax.add_collection(lc)

        # Add scatter points with gradient colors
        scatter = ax.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5, marker='o')

    # Auto-scale to fit LineCollection
    ax.autoscale()
    
    # Set labels for each subplot
    ax.set_xlabel("Green light intensity %", color="green")
    ax.set_ylabel("GFP 395nm")
    ax.set_title(f"Position {i+1}")
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    

# Hide any unused subplots
for j in range(n_positions, len(axs)):
    axs[j].axis('off')

fig.suptitle("GFP 395nm by intensity - all positions", fontsize=14, y=1.00)
fig.tight_layout()
plt.show()



#GFP/OD600 by intensity - all
by_intensity_df = sorted_data_df[["green_intensity","GFP/OD600_raw"]].copy()
by_intensity_df["GFP/OD600_raw"] = by_intensity_df["GFP/OD600_raw"].apply(lambda x: [array[-1] for array in x])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jbl_new_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137new')]
jbl_hybrid_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137hybrid')]
jbl_kirill_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137kirill')]
jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
media_data = by_intensity_df[by_intensity_df.index.str.contains('media')]

colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

fig, axs = plt.subplots()

for i in range(len(jbl_kirill_data["GFP/OD600_raw"].iloc[0])):

    # JCCO plot with gradient color
    x = jbl_kirill_data["green_intensity_percentage"].values
    y = np.array([row[i] for row in jbl_kirill_data['GFP/OD600_raw']])

    # Create line segments for gradient
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Normalize x values for colormap
    norm = plt.Normalize(x.min(), x.max())
    lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
    lc.set_array(x)
    axs.add_collection(lc)

    # Add scatter points with gradient colors
    label_jcco = f'JCCO' if i == 0 else None   
    scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5, marker='o')

""" for i in range(len(jbl_data["GFP/OD600_raw"].iloc[0])):
    axs.plot(jbl_data["green_intensity_percentage"], [row[i] for row in jbl_data['GFP/OD600_raw']],
                label = "JBL001" if i == 0 else None,
                color = "black",
                marker = "^", markersize = 3.0,
                linestyle = "solid", linewidth = 1.0,
                ) """

x_axis = axs.set_xlabel("Green light intensity %")
x_axis.set_color("green")

axs.set_ylabel("GFP/OD600")
axs.set_title(f"GFP/OD600 by intensity - all")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()


#GFP by intensity - all positions in separate subplots

by_intensity_df = sorted_data_df[["green_intensity","GFP/OD600_raw"]].copy()
by_intensity_df["GFP/OD600_raw"] = by_intensity_df["GFP/OD600_raw"].apply(lambda x: [array[-1] for array in x])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jbl_new_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137new')]
jbl_hybrid_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137hybrid')]
jbl_kirill_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137kirill')]
jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
media_data = by_intensity_df[by_intensity_df.index.str.contains('media')]

colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

# Get max number of positions across all datasets
n_jcco = len(jbl_kirill_data["GFP/OD600_raw"].iloc[0]) if len(jbl_kirill_data) > 0 else 0
#n_jbl = len(jbl_data["GFP/OD600_raw"].iloc[0]) if len(jbl_data) > 0 else 0
#n_media = len(media_data["GFP/OD600_raw"].iloc[0]) if len(media_data) > 0 else 0
n_positions = n_jcco

# Create subplots - adjust rows and cols as needed
n_cols = 2  # 2 columns
n_rows = (n_positions + n_cols - 1) // n_cols  # Calculate rows needed
fig, axs = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))
axs = axs.flatten()  # Flatten to 1D array for easier indexing

for i in range(n_positions):
    ax = axs[i]
    
    # JCCO plot with gradient color
    if i < n_jcco:
        x = jbl_kirill_data["green_intensity_percentage"].values
        y = np.array([row[i] for row in jbl_kirill_data['GFP/OD600_raw']])

        # Create line segments for gradient
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Normalize x values for colormap
        norm = plt.Normalize(x.min(), x.max())
        lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
        lc.set_array(x)
        ax.add_collection(lc)

        # Add scatter points with gradient colors
        scatter = ax.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5, marker='o')

    # Auto-scale to fit LineCollection
    ax.autoscale()
    
    # Set labels for each subplot
    ax.set_xlabel("Green light intensity %", color="green")
    ax.set_ylabel("GFP/OD600")
    ax.set_title(f"Position {i+1}")
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    

# Hide any unused subplots
for j in range(n_positions, len(axs)):
    axs[j].axis('off')

fig.suptitle("GFP/OD600 by intensity - all positions", fontsize=14, y=1.00)
fig.tight_layout()
plt.show()



####
#GFP by intensity - average over time

by_intensity_df = sorted_data_df[["green_intensity","GFP/OD600_average","GFP/OD600_std"]].copy()
by_intensity_df["GFP/OD600_average_t0"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[0])
by_intensity_df["GFP/OD600_average_t4"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[1])
by_intensity_df["GFP/OD600_average_t12"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[2])
by_intensity_df["GFP/OD600_average_t24"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[3])
by_intensity_df["GFP/OD600_std_t0"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[0])
by_intensity_df["GFP/OD600_std_t4"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[1])
by_intensity_df["GFP/OD600_std_t12"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[2])
by_intensity_df["GFP/OD600_std_t24"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[3])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jbl_new_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137new')]
jbl_hybrid_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137hybrid')]
jbl_kirill_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137kirill')]
jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
media_data = by_intensity_df[by_intensity_df.index.str.contains('media')]

colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

fig, axs = plt.subplots()

err_color = {"t0": "red",
             "t4": "orange",
             "t12": "green",
             "t24": "blue"}

for i in ["t0","t4","t12","t24"]:

    axs.errorbar(jbl_kirill_data["green_intensity_percentage"],
                 jbl_kirill_data["GFP/OD600_average"+"_"+i],
                yerr = jbl_kirill_data["GFP/OD600_std"+"_"+i], capsize = 2.0,
                color = "black",
                ecolor=err_color[i],
                marker = "^", markersize = 3.0,
                linestyle = "solid", linewidth = 1.0,
                )

    # JCCO plot with gradient color
    x = jbl_kirill_data["green_intensity_percentage"].values
    y = jbl_kirill_data["GFP/OD600_average"+"_"+i].values
    yerr = jbl_kirill_data["GFP/OD600_std"+"_"+i].values

    # Create line segments for gradient
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Normalize x values for colormap
    norm = plt.Normalize(x.min(), x.max())
    lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
    lc.set_array(x)
    axs.add_collection(lc)

    # Add scatter points with gradient colors
    scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5, marker='o')

    # Add error bars for JCCO
    axs.errorbar(x, y, yerr=yerr, fmt='none', ecolor=err_color[i], alpha=1.0, capsize=2.0, label = i)

x_axis = axs.set_xlabel("Green light intensity %")
x_axis.set_color("green")

axs.set_ylabel("GFP 395nm/OD600")
axs.set_title(f"GFP 395nm/OD600 by intensity - average over time")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()


####
#GFP/OD by intensity - average over cell type at 24hr

by_intensity_df = sorted_data_df[["green_intensity","GFP/OD600_average","GFP/OD600_std"]].copy()
by_intensity_df["GFP/OD600_average"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[-1])
by_intensity_df["GFP/OD600_std"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[-1])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jbl_new_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137new')]
jbl_hybrid_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137hybrid')]
jbl_kirill_data = by_intensity_df[by_intensity_df.index.str.contains('JBL137kirill')]
jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
media_data = by_intensity_df[by_intensity_df.index.str.contains('media')]

colors = [(1, 0, 0), (0, 1, 0)]  # Red (1,0,0) to Green (0,1,0)
cmap = LinearSegmentedColormap.from_list('red_green', colors, N=256)

fig, axs = plt.subplots()

err_color = {"JBL137_new": "red",
             "JBL137_hybrid": "orange",
             "JBL137_kirills": "green",
             "JCCO": "blue"}

cell_df = {"JBL137_new": jbl_new_data,
             "JBL137_hybrid": jbl_hybrid_data,
             "JBL137_kirills": jbl_kirill_data,
             "JCCO": jcco_data}

for i in ["JBL137_new","JBL137_hybrid","JBL137_kirills","JCCO"]:

    axs.errorbar(cell_df[i]["green_intensity_percentage"],
                 cell_df[i]["GFP/OD600_average"],
                yerr = cell_df[i]["GFP/OD600_std"], capsize = 2.0,
                color = err_color[i],
                ecolor= err_color[i],
                marker = "^", markersize = 3.0,
                linestyle = "solid", linewidth = 1.0,
                label = i,
                )

    # JCCO plot with gradient color
    x = cell_df[i]["green_intensity_percentage"].values
    y = cell_df[i]["GFP/OD600_average"].values
    yerr = cell_df[i]["GFP/OD600_std"].values

    # Create line segments for gradient
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Normalize x values for colormap
    norm = plt.Normalize(x.min(), x.max())
    lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=1.0)
    lc.set_array(x)
    axs.add_collection(lc)

    # Add scatter points with gradient colors
    scatter = axs.scatter(x, y, c=x, cmap=cmap, norm=norm, s=20, zorder=5, marker='o')

    # Add error bars for JCCO
    axs.errorbar(x, y, yerr=yerr, fmt='none', ecolor=err_color[i], alpha=1.0, capsize=2.0, label = i)

x_axis = axs.set_xlabel("Green light intensity %")
x_axis.set_color("green")

axs.set_ylabel("GFP 395nm/OD600")
axs.set_title(f"GFP 395nm/OD600 by intensity - average over cell type")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()