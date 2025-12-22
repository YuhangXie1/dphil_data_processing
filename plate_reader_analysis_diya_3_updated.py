import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import copy
import json
from datetime import datetime
from os import listdir, path
import re
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection

#functions
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


### MAIN ###
key_rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
key_columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
key_wells = [str(row)+str(col) for row in key_rows for col in key_columns]
plate_map = {key : np.nan for key in key_wells}


NCOLS = 12
NROWS = 8
DEFAULT_OPTICAL_POWER = np.array([

	# Channel 0 (Color 0 or 4). Blue on v0.4c
	[[0 / (2**(7-row)) for col in range(NCOLS)] for row in range(NROWS)],

    # Channel 1 – Green light (checkerboard starting green at B2)
    np.array([
        # Columns 1–12 (A–L), Rows A–H
        [0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # Row A
        [0.0,   2.8,    0.0,    2.0,    0.0,    1.5,    0.0,    1.0,    0.0,    0.5,    0.0,    0.0],  # Row B
        [0.0,   0.0,    1.0,    0.0,    0.5,    0.0,    2.8,    0.0,    2.0,    0.0,    1.5,    0.0],  # Row C
        [0.0,   1.5,    0.0,    2.0,    0.0,    2.8,    0.0,    0.5,    0.0,    1.0,    0.0,    0.0],  # Row D
        [0.0,   0.0,    0.5,    0.0,    1.0,    0.0,    1.5,    0.0,    2.0,    0.0,    2.8,    0.0],  # Row E
        [0.0,   2.0,    0.0,    2.8,    0.0,    0.5,    0.0,    1.0,    0.0,    1.5,    0.0,    0.0],  # Row F
        [0.0,   0.0,    1.0,    0.0,    1.5,    0.0,    2.0,    0.0,    2.8,    0.0,    0.5,    0.0],  # Row Gh
        [0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # Row H
    ]),
	# Channel 2 (Color 2 or 6). Yellow-Green or White on v0.4c
	[[0 / (2**row) for col in range(NCOLS)] for row in range(NROWS)],    

    # Channel 3 – Red light (opposite checkerboard cells)
    np.array([
        [0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # Row A
        [0.0,   0.0,    1.0,    0.0,    0.5,    0.0,    2.8,    0.0,    2.0,    0.0,    1.5,    0.0],  # Row B
        [0.0,   2.8,    0.0,    2.0,    0.0,    1.5,    0.0,    1.0,    0.0,    0.5,    0.0,    0.0],  # Row C
        [0.0,   0.0,    0.5,    0.0,    1.0,    0.0,    1.5,    0.0,    2.0,    0.0,    2.8,    0.0],  # Row D
        [0.0,   1.5,    0.0,    2.0,    0.0,    2.8,    0.0,    0.5,    0.0,    1.0,    0.0,    0.0],  # Row E
        [0.0,   0.0,    1.0,    0.0,    1.5,    0.0,    2.0,    0.0,    2.8,    0.0,    0.5,    0.0],  # Row F
        [0.0,   2.0,    0.0,    2.8,    0.0,    0.5,    0.0,    1.0,    0.0,    1.5,    0.0,    0.0],  # Row G
        [0.0,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0],  # Row H
    ])
])

#turn light array into labels - need to code
plate_map_new = {key : np.nan for key in key_wells}
for key in plate_map_new:
    for channel in DEFAULT_OPTICAL_POWER:
        
        pass



#plate map
plate_map = {
    # Row A
    "A1": "media_dark_-1",
    "A2": "JBL001_dark_-1", "A3": "JBL001_dark_-1", "A4": "JBL001_dark_-1",
    "A5": "JBL001_dark_-1", "A6": "JBL001_dark_-1", "A7": "JBL001_dark_-1", "A8": "JBL001_dark_-1",
    "A9": "JBL001_dark_-1", "A10": "JBL001_dark_-1", "A11": "JBL001_dark_-1",
    "A12": "media_dark_-1",

    # Row B 
    "B1": "JBL001_dark_-1",
    "B2": "JBL36_green_2.8",  "B3": "JBL36_red_1.0",  "B4": "JBL36_green_2.0",
    "B5": "JBL36_red_0.5",    "B6": "JBL36_green_1.5", "B7": "JBL36_red_2.8",
    "B8": "JBL36_green_1.0",  "B9": "JBL36_red_2.0",  "B10": "JBL36_green_0.5",
    "B11": "JBL36_red_1.5",  
    "B12": "JBL001_dark_-1",

    # Row C
    "C1": "JBL001_dark_-1",
    "C2": "JCCO_red_2.8",   "C3": "JCCO_green_1.0",  "C4": "JCCO_red_2.0",
    "C5": "JCCO_green_0.5", "C6": "JCCO_red_1.5",    "C7": "JCCO_green_2.8",
    "C8": "JCCO_red_1.0",   "C9": "JCCO_green_2.0",  "C10": "JCCO_red_0.5",
    "C11": "JCCO_green_1.5",
    "C12": "JBL001_dark_-1",

    # Row D
    "D1": "JBL001_dark_-1",
    "D2": "JBL36_green_1.5",  "D3": "JBL36_red_0.5",  "D4": "JBL36_green_2.0",
    "D5": "JBL36_red_1.0",    "D6": "JBL36_green_2.8", "D7": "JBL36_red_1.5",
    "D8": "JBL36_green_0.5",  "D9": "JBL36_red_2.0",  "D10": "JBL36_green_1.0",
    "D11": "JBL36_red_2.8",  
    "D12": "JBL001_dark_-1",

    # Row E
    "E1": "JBL001_dark_-1",
    "E2": "JCCO_red_1.5",    "E3": "JCCO_green_0.5",  "E4": "JCCO_red_2.0",
    "E5": "JCCO_green_1.0",  "E6": "JCCO_red_2.8",    "E7": "JCCO_green_1.5",
    "E8": "JCCO_red_0.5",    "E9": "JCCO_green_2.0",  "E10": "JCCO_red_1.0",
    "E11": "JCCO_green_2.8",
    "E12": "JBL001_dark_-1",

    # Row F
    "F1": "JBL001_dark_-1",
    "F2": "JBL36_green_2.0",  "F3": "JBL36_red_1.0",  "F4": "JBL36_green_2.8",
    "F5": "JBL36_red_1.5",    "F6": "JBL36_green_0.5", "F7": "JBL36_red_2.0",
    "F8": "JBL36_green_1.0",  "F9": "JBL36_red_2.8",  "F10": "JBL36_green_1.5",
    "F11": "JBL36_red_0.5",  
    "F12": "JBL001_dark_-1",

    # Row G 
    "G1": "JBL001_dark_-1",
    "G2": "JCCO_red_2.0",    "G3": "JCCO_green_1.0",  "G4": "JCCO_red_2.8",
    "G5": "JCCO_green_1.5",  "G6": "JCCO_red_0.5",    "G7": "JCCO_green_2.0",
    "G8": "JCCO_red_1.0",    "G9": "JCCO_green_2.8",  "G10": "JCCO_red_1.5",
    "G11": "JCCO_green_0.5", 
    "G12": "JBL001_dark_-1",

    # Row H
    "H1": "media_dark_-1",
    "H2": "JBL001_dark_-1", "H3": "JBL001_dark_-1", "H4": "JBL001_dark_-1",
    "H5": "JBL001_dark_-1", "H6": "JBL001_dark_-1", "H7": "JBL001_dark_-1", "H8": "JBL001_dark_-1",
    "H9": "JBL001_dark_-1", "H10": "JBL001_dark_-1", "H11": "JBL001_dark_-1",
    "H12": "media_dark_-1",

}
#Note JCCO_green_3.2 is actually green 2.8 with no red light. Rest have 2.8 red light

#color and linestyles
#color_map = {"green": "green",
#             "red": "red",
#}
color_map = {"2.8": (0.0,1.0,0.0),
                 "1.4": (0.1,0.9,0.0),
                 "0.84": (0.2,0.8,0.0),
                 "0.56": (0.3,0.7,0.0),
                 "0.28": (0.4,0.6,0.0),
                 "0.14": (0.5,0.5,0.0),
                 "0.07":(0.6,0.4,0.0),
                 "0.028":(0.7,0.3,0.0),
                 "0.014":(0.8,0.2,0.0),
                 "0.0028":(0.9,0.1,0.0),
                 "0":(1.0,0.0,0.0),
                 "-1":"black",
                 "":"black",
                 "3.2":"darkgreen",
}

color_map = {"red":"red",
             "green":"green",
             "dark":"black",
}

linestyle_map = {"JBL137": "solid",
                 "JBL36": "solid",
                 "JBL001":"dashed",
                 "JCCO":"solid",
                 "media":"dotted",
}

alpha_map = {"2.8": 1,
                 "2.0": 0.9,
                 "1.5": 0.8,
                 "1.0": 0.6,
                 "0.5": 0.4,
                 "0": 0.2,
                 "-1": 1,
                 "": 1,
                 "3.2": 1,
}

marker_map = {"JBL137": "o",
              "JBL36": "^",
              "JBL001":"^",
              "JCCO":"o",
              "media":"^",
}

#filepaths = ["25-11-19_diya_4/25-11-18_diya2_dose_curve_low_gain_extracted_OD600.csv","25-11-19_diya_4/25-11-18_diya2_dose_curve_low_gain_extracted_GFP 475nm.csv", "25-11-19_diya_4/25-11-18_diya2_dose_curve_low_gain_extracted_GFP 395nm.csv"]

filepaths = ["25-11-06_diya_3/Diya/25-11-06_diya_3_diya_extracted_OD600.csv", "25-11-06_diya_3/Diya/25-11-06_diya_3_diya_extracted_GFP 395nm.csv"]


#initiating data df
indexes = sorted(set(plate_map.values()))
sorted_data_df = pd.DataFrame(0.0, index=indexes, columns=["green_intensity",
                                                           "red_intensity",
                                                           "timepoints",
                                                           "OD600_raw_array","OD600_average","OD600_std",
                                                           #"GFP475_raw_array","GFP475_average","GFP475_std",
                                                           "GFP395_raw_array","GFP395_average","GFP395_std",
                                                           "GFP/OD600_raw","GFP/OD600_average","GFP/OD600_std"]
                                                            ).astype(object)
#initialising arrays
sorted_data_df["OD600_raw_array"] = [[] for _ in range(len(sorted_data_df))]
#sorted_data_df["GFP475_raw_array"] = [[] for _ in range(len(sorted_data_df))]
sorted_data_df["GFP395_raw_array"] = [[] for _ in range(len(sorted_data_df))]
sorted_data_df["GFP/OD600_raw"] = [[] for _ in range(len(sorted_data_df))]

file_list = {}
for filepath in filepaths:
    filename = filepath.split("_")[-1].replace(".csv","")
    data = pd.read_csv(filepath, header = 0, index_col= 0)
    
    #working out the time points
    timepoints = list(data.columns.values)
    try:
        initial_time = datetime.strptime(timepoints[0], "%Y-%m-%d %H:%M:%S")
    except ValueError:
        try:
            initial_time = datetime.strptime(timepoints[0], "%d/%m/%Y %H:%M:%S")
        except ValueError:
            initial_time = datetime.strptime(timepoints[0], "%m/%d/%Y %I:%M:%S %p")

    timepoints_hrs = []
    for item in timepoints:
        try:
            timepoints_datetime = datetime.strptime(item, "%Y-%m-%d %H:%M:%S")
        except ValueError:
            try:
                timepoints_datetime = datetime.strptime(item, "%d/%m/%Y %H:%M:%S")
            except ValueError:
                timepoints_datetime = datetime.strptime(item, "%m/%d/%Y %I:%M:%S %p")
        time_delta = timepoints_datetime - initial_time
        timepoints_hrs.append(time_delta.total_seconds() / 3600)

    #renaming file column to the time
    data = data.rename(columns={old:new for (old,new) in zip(timepoints, timepoints_hrs)})
    
    file_list[filename] = data


#populate raw data from wells
for well_coord, name in plate_map.items():
    sorted_data_df.loc[name,"OD600_raw_array"].append(np.array(file_list["OD600"].loc[str(well_coord)]))
    #sorted_data_df.loc[name,"GFP475_raw_array"].append(np.array(file_list["GFP 475nm"].loc[str(well_coord)]))
    sorted_data_df.loc[name,"GFP395_raw_array"].append(np.array(file_list["GFP 395nm"].loc[str(well_coord)]))
    sorted_data_df.loc[name,"GFP/OD600_raw"].append(np.array(file_list["GFP 395nm"].loc[str(well_coord)])/np.array(file_list["OD600"].loc[str(well_coord)]))

#intensity and time points
sorted_data_df["timepoints"] = [list(file_list["OD600"].columns.values) for _ in range(len(sorted_data_df))]
for index, row in sorted_data_df.iterrows():
    sorted_data_df.loc[index,"green_intensity"] = float(index.split("_")[-1])
    sorted_data_df.loc[index,"red_intensity"] = float(index.split("_")[-1])

    #mean, std data
    sorted_data_df.at[index,"OD600_average"] = np.mean(sorted_data_df.loc[index,"OD600_raw_array"], axis = 0)
    sorted_data_df.at[index,"OD600_std"] = np.std(sorted_data_df.loc[index,"OD600_raw_array"], axis = 0)
    #sorted_data_df.at[index,"GFP475_average"] = np.mean(sorted_data_df.loc[index,"GFP475_raw_array"], axis = 0)
    #sorted_data_df.at[index,"GFP475_std"] = np.std(sorted_data_df.loc[index,"GFP475_raw_array"], axis = 0)
    sorted_data_df.at[index,"GFP395_average"] = np.mean(sorted_data_df.loc[index,"GFP395_raw_array"], axis = 0)
    sorted_data_df.at[index,"GFP395_std"] = np.std(sorted_data_df.loc[index,"GFP395_raw_array"], axis = 0)

    sorted_data_df.at[index,"GFP/OD600_average"] = np.mean(sorted_data_df.loc[index,"GFP/OD600_raw"], axis = 0)
    sorted_data_df.at[index,"GFP/OD600_std"] = np.std(sorted_data_df.loc[index,"GFP/OD600_raw"], axis = 0)







#plotting

plot_selection = {
    "cells":["JBL137","JBL36"],
    "color":["green","red"],
    "intensity":[2.8,2.3,1.8,1.3,0.8,0.3,0]
}

plot_exclude = {
    "cells":["JBL36","media"],
    "color":[],
    "intensity":[]
}

#od600 - average
fig, axs = plt.subplots()
for row_index, row in sorted_data_df.iterrows():
    index_name_cells, index_name_color, index_name_intensity = row_index.split("_")
    if index_name_cells in plot_exclude["cells"] or index_name_color in plot_exclude["color"] or index_name_intensity in plot_exclude["intensity"]:
        continue
    else:
        axs.errorbar(row["timepoints"], row["OD600_average"],
                    yerr = row["OD600_std"], capsize = 2.0,
                    label = row_index,
                    color = color_map[index_name_color],
                    marker = marker_map[index_name_cells], markersize = 3.0,
                    linestyle = linestyle_map[index_name_cells], linewidth = 1.0,
                    alpha = alpha_map[index_name_intensity]
                    )
#alpha = alpha_map[index_name_intensity]
axs.set_xlabel("Time (hrs)")
axs.set_ylabel("OD600")
axs.set_title(f"OD600 - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#GFP - average
fig, axs = plt.subplots()
for row_index, row in sorted_data_df.iterrows():
    index_name_cells, index_name_color, index_name_intensity = row_index.split("_")
    if index_name_cells in plot_exclude["cells"] or index_name_color in plot_exclude["color"] or index_name_intensity in plot_exclude["intensity"]:
        continue
    else:
        axs.errorbar(row["timepoints"], row["GFP395_average"],
                    yerr = row["GFP395_std"], capsize = 2.0,
                    label = row_index,
                    color = color_map[index_name_color],
                    marker = marker_map[index_name_cells], markersize = 3.0,
                    linestyle = linestyle_map[index_name_cells], linewidth = 1.0,
                    alpha = alpha_map[index_name_intensity]
                    )
#alpha = alpha_map[index_name_intensity]
axs.set_xlabel("Time (hrs)")
axs.set_ylabel("GFP 395nm")
axs.set_title(f"GFP 395nm - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#GFP - average
""" fig, axs = plt.subplots()
for row_index, row in sorted_data_df.iterrows():
    index_name_cells, index_name_color, index_name_intensity = row_index.split("_")
    if index_name_cells in plot_exclude["cells"] or index_name_color in plot_exclude["color"] or index_name_intensity in plot_exclude["intensity"]:
        continue
    else:
        axs.errorbar(row["timepoints"], row["GFP475_average"],
                    yerr = row["GFP475_std"], capsize = 2.0,
                    label = row_index,
                    color = color_map[index_name_intensity],
                    marker = marker_map[index_name_cells], markersize = 3.0,
                    linestyle = linestyle_map[index_name_cells], linewidth = 1.0,
                    )
#alpha = alpha_map[index_name_intensity]
axs.set_xlabel("Time (hrs)")
axs.set_ylabel("GFP 475nm")
axs.set_title(f"GFP 475nm - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show() """

#gfp/od600 - average
fig, axs = plt.subplots()
for row_index, row in sorted_data_df.iterrows():
    index_name_cells, index_name_color, index_name_intensity = row_index.split("_")
    if index_name_cells in plot_exclude["cells"] or index_name_color in plot_exclude["color"] or index_name_intensity in plot_exclude["intensity"]:
        continue
    else:
        axs.errorbar(row["timepoints"], row["GFP/OD600_average"],
                    yerr = row["GFP/OD600_std"], capsize = 2.0,
                    label = row_index,
                    color = color_map[index_name_color],
                    marker = marker_map[index_name_cells], markersize = 3.0,
                    linestyle = linestyle_map[index_name_cells], linewidth = 1.0,
                    alpha = alpha_map[index_name_intensity],
                    )
#alpha = alpha_map[index_name_intensity]
axs.set_xlabel("Time (hrs)")
axs.set_ylabel("GFP 395/OD600")
axs.set_title(f"GFP 395/OD600 - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()


#OD600 - all lines
fig, axs = plt.subplots()
for row_index, row in sorted_data_df.iterrows():
    index_name_cells, index_name_color, index_name_intensity = row_index.split("_")
    if index_name_cells in plot_exclude["cells"] or index_name_color in plot_exclude["color"] or index_name_intensity in plot_exclude["intensity"]:
        continue
    if index_name_intensity in []:
        #["2.8","2.52","2.24","1.96","1.68","1.4","1.12","0.84","0.56","0.28","0","-1"]
        continue
    else:
        for repeat in row["OD600_raw_array"]:
            axs.plot(row["timepoints"], repeat,
                        label = row_index,
                        color = color_map[index_name_color],
                        marker = marker_map[index_name_cells], markersize = 3.0,
                        linestyle = linestyle_map[index_name_cells], linewidth = 1.0,
                        alpha = alpha_map[index_name_intensity],
                        )
#alpha = alpha_map[index_name_intensity]
axs.set_xlabel("Time (hrs)")
axs.set_ylabel("OD600")
axs.set_title(f"OD600 - all")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#GFP 395 - all lines
fig, axs = plt.subplots()
for row_index, row in sorted_data_df.iterrows():
    index_name_cells, index_name_color, index_name_intensity = row_index.split("_")
    if index_name_cells in plot_exclude["cells"] or index_name_color in plot_exclude["color"] or index_name_intensity in plot_exclude["intensity"]:
        continue
    if index_name_intensity in []:
        #["2.8","2.52","2.24","1.96","1.68","1.4","1.12","0.84","0.56","0.28","0","-1"]
        continue
    else:
        for repeat in row["GFP395_raw_array"]:
            axs.plot(row["timepoints"], repeat,
                        label = row_index,
                        color = color_map[index_name_color],
                        marker = marker_map[index_name_cells], markersize = 3.0,
                        linestyle = linestyle_map[index_name_cells], linewidth = 1.0,
                        alpha = alpha_map[index_name_intensity],
                        )
#alpha = alpha_map[index_name_intensity]
axs.set_xlabel("Time (hrs)")
axs.set_ylabel("GFP 395nm")
axs.set_title(f"GFP 395nm - all")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#GFP 395/OD600 - all lines
fig, axs = plt.subplots()
for row_index, row in sorted_data_df.iterrows():
    index_name_cells, index_name_color, index_name_intensity = row_index.split("_")
    if index_name_cells in plot_exclude["cells"] or index_name_color in plot_exclude["color"] or index_name_intensity in plot_exclude["intensity"]:
        continue
    if index_name_intensity in []:
        #["2.8","2.52","2.24","1.96","1.68","1.4","1.12","0.84","0.56","0.28","0","-1"]
        continue
    else:
        for repeat in row["GFP/OD600_raw"]:
            axs.plot(row["timepoints"], repeat,
                        label = row_index,
                        color = color_map[index_name_color],
                        marker = marker_map[index_name_cells], markersize = 3.0,
                        linestyle = linestyle_map[index_name_cells], linewidth = 1.0,
                        alpha = alpha_map[index_name_intensity],
                        )
#alpha = alpha_map[index_name_intensity]
axs.set_xlabel("Time (hrs)")
axs.set_ylabel("GFP 395nm/OD600")
axs.set_title(f"GFP 395nm/OD600 - all")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()



#GFP by intensity - average

by_intensity_df = sorted_data_df[["green_intensity","GFP395_average","GFP395_std"]].copy()
by_intensity_df["GFP395_average"] = by_intensity_df["GFP395_average"].apply(lambda x: x[-2])
by_intensity_df["GFP395_std"] = by_intensity_df["GFP395_std"].apply(lambda x: x[-2])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
jcco_green_data = jcco_data[jcco_data.index.str.contains('green')]
jcco_red_data = jcco_data[jcco_data.index.str.contains('red')]

fig, axs = plt.subplots()

axs.errorbar(jcco_green_data["green_intensity_percentage"], jcco_green_data["GFP395_average"],
            yerr = jcco_green_data["GFP395_std"], capsize = 2.0,
            label = "JCCO - green",
            color = "green",
            marker = "o", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            )

axs.errorbar(jcco_red_data["green_intensity_percentage"], jcco_red_data["GFP395_average"],
            yerr = jcco_red_data["GFP395_std"], capsize = 2.0,
            label = "JCCO - red",
            color = "red",
            marker = "o", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            )


axs.set_xlabel("Light intensity %")
axs.set_ylabel("GFP 395nm")
axs.set_title(f"GFP 395nm by intensity - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#OD600 by intensity - average

by_intensity_df = sorted_data_df[["green_intensity","OD600_average","OD600_std"]].copy()
by_intensity_df["OD600_average"] = by_intensity_df["OD600_average"].apply(lambda x: x[-2])
by_intensity_df["OD600_std"] = by_intensity_df["OD600_std"].apply(lambda x: x[-2])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
jcco_green_data = jcco_data[jcco_data.index.str.contains('green')]
jcco_red_data = jcco_data[jcco_data.index.str.contains('red')]

fig, axs = plt.subplots()

axs.errorbar(jcco_green_data["green_intensity_percentage"], jcco_green_data["OD600_average"],
            yerr = jcco_green_data["OD600_std"], capsize = 2.0,
            label = "JCCO - green",
            color = "green",
            marker = "o", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            )

axs.errorbar(jcco_red_data["green_intensity_percentage"], jcco_red_data["OD600_average"],
            yerr = jcco_red_data["OD600_std"], capsize = 2.0,
            label = "JCCO - red",
            color = "red",
            marker = "o", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            )


axs.set_xlabel("Light intensity %")
axs.set_ylabel("OD600")
axs.set_title(f"OD600 by intensity - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#GFP/OD600 by intensity - average

by_intensity_df = sorted_data_df[["green_intensity","GFP/OD600_average","GFP/OD600_std"]].copy()
by_intensity_df["GFP/OD600_average"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[-1])
by_intensity_df["GFP/OD600_std"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[-1])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
jcco_green_data = jcco_data[jcco_data.index.str.contains('green')]
jcco_red_data = jcco_data[jcco_data.index.str.contains('red')]

fig, axs = plt.subplots()

axs.errorbar(jcco_green_data["green_intensity_percentage"], jcco_green_data["GFP/OD600_average"],
            yerr = jcco_green_data["GFP/OD600_std"], capsize = 2.0,
            label = "JCCO - green",
            color = "green",
            marker = "o", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            )

axs.errorbar(jcco_red_data["green_intensity_percentage"], jcco_red_data["GFP/OD600_average"],
            yerr = jcco_red_data["GFP/OD600_std"], capsize = 2.0,
            label = "JCCO - red",
            color = "red",
            marker = "o", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            )


axs.set_xlabel("Light intensity %")
axs.set_ylabel("GFP 395nm/OD600")
axs.set_title(f"GFP 395nm/OD600 by intensity - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()



#GFP by intensity - all
by_intensity_df = sorted_data_df[["green_intensity","GFP395_raw_array"]].copy()
by_intensity_df["GFP395_raw_array"] = by_intensity_df["GFP395_raw_array"].apply(lambda x: [array[-1] for array in x])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
jcco_green_data = jcco_data[jcco_data.index.str.contains('green')]
jcco_red_data = jcco_data[jcco_data.index.str.contains('red')]

fig, axs = plt.subplots()

for i in range(len(jcco_data["GFP395_raw_array"].iloc[0])):

    axs.plot(jcco_green_data["green_intensity_percentage"], np.array([row[i] for row in jcco_green_data['GFP395_raw_array']]),
                label = "JCCO - green",
                color = "green",
                marker = "o", markersize = 3.0,
                linestyle = "solid", linewidth = 1.0,
                )
    
    axs.plot(jcco_red_data["green_intensity_percentage"], np.array([row[i] for row in jcco_red_data['GFP395_raw_array']]),
            label = "JCCO - red",
            color = "red",
            marker = "o", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            )

axs.set_xlabel("Light intensity %")
axs.set_ylabel("GFP 395nm")
axs.set_title("GFP 395nm by intensity - all")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#GFP/OD600 by intensity - all
by_intensity_df = sorted_data_df[["green_intensity","GFP/OD600_raw"]].copy()
by_intensity_df["GFP/OD600_raw"] = by_intensity_df["GFP/OD600_raw"].apply(lambda x: [array[-1] for array in x])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
jcco_green_data = jcco_data[jcco_data.index.str.contains('green')]
jcco_red_data = jcco_data[jcco_data.index.str.contains('red')]

fig, axs = plt.subplots()

for i in range(len(jcco_data["GFP/OD600_raw"].iloc[0])):

    axs.plot(jcco_green_data["green_intensity_percentage"], np.array([row[i] for row in jcco_green_data["GFP/OD600_raw"]]),
                label = "JCCO - green",
                color = "green",
                marker = "o", markersize = 3.0,
                linestyle = "solid", linewidth = 1.0,
                )
    
    axs.plot(jcco_red_data["green_intensity_percentage"], np.array([row[i] for row in jcco_red_data["GFP/OD600_raw"]]),
            label = "JCCO - red",
            color = "red",
            marker = "o", markersize = 3.0,
            linestyle = "solid", linewidth = 1.0,
            )

axs.set_xlabel("Light intensity %")
axs.set_ylabel("GFP 395nm/OD600")
axs.set_title("GFP 395nm/OD600 by intensity - all")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()



####
#GFP/OD by intensity - average

by_intensity_df = sorted_data_df[["green_intensity","GFP/OD600_average","GFP/OD600_std"]].copy()
by_intensity_df["GFP/OD600_average_t0"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[0])
by_intensity_df["GFP/OD600_average_t4"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[1])
by_intensity_df["GFP/OD600_average_t12"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[2])
by_intensity_df["GFP/OD600_average_t27"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[3])
by_intensity_df["GFP/OD600_average_t49"] = by_intensity_df["GFP/OD600_average"].apply(lambda x: x[4])
by_intensity_df["GFP/OD600_std_t0"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[0])
by_intensity_df["GFP/OD600_std_t4"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[1])
by_intensity_df["GFP/OD600_std_t12"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[2])
by_intensity_df["GFP/OD600_std_t27"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[3])
by_intensity_df["GFP/OD600_std_t49"] = by_intensity_df["GFP/OD600_std"].apply(lambda x: x[4])
by_intensity_df["green_intensity_percentage"] = by_intensity_df["green_intensity"].apply(lambda x: 100*x/2.8)
by_intensity_df = by_intensity_df[~by_intensity_df.index.str.contains('-1')]

jcco_data = by_intensity_df[by_intensity_df.index.str.contains('JCCO')]
jcco_green_data = jcco_data[jcco_data.index.str.contains('green')]
jcco_red_data = jcco_data[jcco_data.index.str.contains('red')]


fig, axs = plt.subplots()

err_color = {"t0": "red",
             "t4": "orange",
             "t12": "green",
             "t27": "blue",
             "t49": "brown"}

for i in err_color.keys():

    axs.errorbar(jcco_green_data["green_intensity_percentage"], jcco_green_data["GFP/OD600_average"+"_"+i],
                yerr = jcco_green_data["GFP/OD600_std"+"_"+i], capsize = 2.0,
                label = "JCCO - green " + i,
                color = "green",
                ecolor=err_color[i],
                marker = "o", markersize = 3.0,
                linestyle = "solid", linewidth = 1.0,
                )
    
    axs.errorbar(jcco_red_data["green_intensity_percentage"], jcco_red_data["GFP/OD600_average"+"_"+i],
                yerr = jcco_red_data["GFP/OD600_std"+"_"+i], capsize = 2.0,
                label = "JCCO - red " + i,
                color = "red",
                ecolor=err_color[i],
                marker = "o", markersize = 3.0,
                linestyle = "solid", linewidth = 1.0,
                )

axs.set_xlabel("Light intensity %")
axs.set_ylabel("GFP 395nm/OD600")
axs.set_title(f"GFP 395nm/OD600 by intensity - average")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()