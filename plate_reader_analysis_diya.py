import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import copy
import json
from datetime import datetime
from os import listdir, path
import re

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
        #parsing name
        namemet_repeat = key.rsplit("_",1)
        if namemet_repeat[0] in data_dict_repeats.keys():
            data_dict_repeats[namemet_repeat[0]].append(value)
        else:
            data_dict_repeats[namemet_repeat[0]] = [value]

    data_dict_processed = {}
    for key, value in data_dict_repeats.items():
        #getting mean, std, cv
        mean_std_cv_array = mean_std_cv(value) #mean_std_cv_array = [[mean_od,mean_gfp, mean_gfp/od],[std_od,std_gfp, std_gfp/od],[cv_od,cv_gfp,cv_gfp/od]]

        #new dictionary with processed data. value = [mean, std, cv]
        data_dict_processed[key] = mean_std_cv_array

    return data_dict_processed

def load_plate_into_df(data, data_type = "xlsx"):
    od600_map = pd.DataFrame(0.0, index = ["A","B","C","D","E","F","G","H"], columns= np.arange(1,13))
    gfp_map = pd.DataFrame(0.0, index = ["A","B","C","D","E","F","G","H"], columns= np.arange(1,13))
    if data_type == "csv":
        for row_num, row in enumerate(data):
            if row_num == 25:
                timepoints_timeformat = row[1]
            if row_num in location_of_data_in_sheet_od600[1]:
                od600_map.loc[row[0]] = np.array(row[1::], dtype=float)
            if row_num in location_of_data_in_sheet_gfp[1]:
                gfp_map.loc[row[0]] = np.array(row[1::], dtype=float)

    elif data_type == "xlsx":
        for row_num, row in enumerate(data.values):
            if row_num == 29:
                timepoints_timeformat = row[1]
            if row_num in location_of_data_in_sheet_od600[1]:
                od600_map.loc[row[0]] = np.array(row[1::], dtype=float)
            if row_num in location_of_data_in_sheet_gfp[1]:
                gfp_map.loc[row[0]] = np.array(row[1::], dtype=float)

    return od600_map, gfp_map, timepoints_timeformat

### MAIN ###

#plate map
plate_map = {
    # Row A - JBL137
    "A1": "JBL137_green_2.8_1", "A2": "JBL137_green_2.3_1", "A3": "JBL137_green_1.8_1", "A4": "JBL137_green_1.3_1",
    "A5": "JBL137_green_0.8_1", "A6": "JBL137_green_0.3_1", "A7": "JBL137_red_0.3_1", "A8": "JBL137_red_0.8_1",
    "A9": "JBL137_red_1.3_1", "A10": "JBL137_red_1.8_1", "A11": "JBL137_red_2.3_1", "A12": "JBL137_red_2.8_1",

    # Row B - JBL36
    "B1": "JBL36_green_2.8_1", "B2": "JBL36_green_2.3_1", "B3": "JBL36_green_1.8_1", "B4": "JBL36_green_1.3_1",
    "B5": "JBL36_green_0.8_1", "B6": "JBL36_green_0.3_1", "B7": "JBL36_red_0.3_1", "B8": "JBL36_red_0.8_1",
    "B9": "JBL36_red_1.3_1", "B10": "JBL36_red_1.8_1", "B11": "JBL36_red_2.3_1", "B12": "JBL36_red_2.8_1",

    # Row C - JBL137
    "C1": "JBL137_green_2.8_2", "C2": "JBL137_green_2.3_2", "C3": "JBL137_green_1.8_2", "C4": "JBL137_green_1.3_2",
    "C5": "JBL137_green_0.8_2", "C6": "JBL137_green_0.3_2", "C7": "JBL137_red_0.3_2", "C8": "JBL137_red_0.8_2",
    "C9": "JBL137_red_1.3_2", "C10": "JBL137_red_1.8_2", "C11": "JBL137_red_2.3_2", "C12": "JBL137_red_2.8_2",

    # Row D - JBL36
    "D1": "JBL36_green_2.8_2", "D2": "JBL36_green_2.3_2", "D3": "JBL36_green_1.8_2", "D4": "JBL36_green_1.3_2",
    "D5": "JBL36_green_0.8_2", "D6": "JBL36_green_0.3_2", "D7": "JBL36_red_0.3_2", "D8": "JBL36_red_0.8_2",
    "D9": "JBL36_red_1.3_2", "D10": "JBL36_red_1.8_2", "D11": "JBL36_red_2.3_2", "D12": "JBL36_red_2.8_2",

    # Row E - JBL137
    "E1": "JBL137_green_2.8_3", "E2": "JBL137_green_2.3_3", "E3": "JBL137_green_1.8_3", "E4": "JBL137_green_1.3_3",
    "E5": "JBL137_green_0.8_3", "E6": "JBL137_green_0.3_3", "E7": "JBL137_red_0.3_3", "E8": "JBL137_red_0.8_3",
    "E9": "JBL137_red_1.3_3", "E10": "JBL137_red_1.8_3", "E11": "JBL137_red_2.3_3", "E12": "JBL137_red_2.8_3",

    # Row F - JBL36
    "F1": "JBL36_green_2.8_3", "F2": "JBL36_green_2.3_3", "F3": "JBL36_green_1.8_3", "F4": "JBL36_green_1.3_3",
    "F5": "JBL36_green_0.8_3", "F6": "JBL36_green_0.3_3", "F7": "JBL36_red_0.3_3", "F8": "JBL36_red_0.8_3",
    "F9": "JBL36_red_1.3_3", "F10": "JBL36_red_1.8_3", "F11": "JBL36_red_2.3_3", "F12": "JBL36_red_2.8_3",

    # Row G - JBL137
    #"G1": "JBL137_green_2.8,1.8_1", "G2": "JBL137_green_2.3,1.8_1", "G3": "JBL137_green_1.8,1.8_1", "G4": "JBL137_green_1.3,1.8_1",
    #"G5": "JBL137_green_0.8,1.8_1", "G6": "JBL137_green_0_1", "G7": "JBL137_red_0_1", "G8": "JBL137_red_1.8,0.8_1",
    #"G9": "JBL137_red_1.8,1.3_1", "G10": "JBL137_red_1.8,1.8_1", "G11": "JBL137_red_1.8,2.3_1", "G12": "JBL137_red_1.8,2.8_1",

    # Row H - JBL36
    #"H1": "JBL36_green_2.8,1.8_1", "H2": "JBL36_green_2.3,1.8_1", "H3": "JBL36_green_1.8,1.8_1", "H4": "JBL36_green_1.3,1.8_1",
    #"H5": "JBL36_green_0.8,1.8_1", "H6": "JBL36_green_0_1", "H7": "JBL36_red_0_1", "H8": "JBL36_red_1.8,0.8_1",
    #"H9": "JBL36_red_1.8,1.3_1", "H10": "JBL36_red_1.8,1.8_1", "H11": "JBL36_red_1.8,2.3_1", "H12": "JBL36_red_1.8,2.8_1"
}

#color and linestyles
color_map = {"green": "green",
             "red": "red",
}

linestyle_map = {"JBL137": "solid",
                 "JBL36": "solid",
}

alpha_map = {"2.8": 1,
                 "2.3": 0.8,
                 "1.8": 0.7,
                 "1.3": 0.5,
                 "0.8": 0.4,
                 "0.3": 0.2,
}

marker_map = {"JBL137": "o",
              "JBL36": "^",
}

#filepaths = ["25-10-07_growth/25-10-07 Dark t0.xlsx","25-10-07_growth/25-10-07 Dark t2.xlsx","25-10-07_growth/25-10-07 Dark t4.xlsx","25-10-07_growth/25-10-07 Dark t19.xlsx","25-10-07_growth/25-10-07 Dark t22.xlsx"]
#filepaths = ["25-10-07_growth/25-10-07 Green t0.xlsx","25-10-07_growth/25-10-07 Green t2.xlsx","25-10-07_growth/25-10-07 Green t4.xlsx","25-10-07_growth/25-10-07 Green t19.xlsx","25-10-07_growth/25-10-07 Green t22.xlsx"]
#filepaths = ["25-10-28_diya/25-10-28_diya_t0.xlsx","25-10-28_diya/25-10-28_diya_t1.xlsx","25-10-28_diya/25-10-28_diya_t1.5.xlsx","25-10-28_diya/25-10-28_diya_t2.xlsx","25-10-28_diya/25-10-28_diya_t3.xlsx"]
filepath = "25-10-28_diya_modified"
filepaths = [path.join(filepath, file) for file in listdir(filepath)]

def extract_time_sorting_key(filename):
    match = re.search(r"_t(\d+(?:\.\d+)?)", filename)
    return float(match.group(1)) if match else -1  # -1 for files without 't'    

filepaths = sorted(filepaths, key=extract_time_sorting_key)


location_of_data_in_sheet_od600 = [np.arange(0,12),np.arange(33,41)] #[columns, rows]
location_of_data_in_sheet_gfp = [np.arange(0,12),np.arange(65,73)] #[columns, rows]

#gathering data from each time point and putting into dict
data_dict = {} #name : [od600, gfp]
timepoints = []
for filepath in filepaths:
    data = pd.read_excel(filepath, 0, header=None, engine="calamine")
    od600_map, gfp_map, timepoints_timeformat = load_plate_into_df(data)
    timepoints.append(timepoints_timeformat)
    
    for key, value in plate_map.items():
        od600_map_item = od600_map.loc[key[0],int(key[1:])]
        gfp_map_item = gfp_map.loc[key[0],int(key[1:])]
        if value not in data_dict.keys():
            data_dict[value] = [[od600_map_item], [gfp_map_item], [gfp_map_item/od600_map_item]]
        else:
            data_dict[value][0].append(od600_map_item)
            data_dict[value][1].append(gfp_map_item)
            data_dict[value][2].append(gfp_map_item/od600_map_item)

#working out the time points
initial_time = datetime.strptime(timepoints[0], "%m/%d/%Y %I:%M:%S %p")
timepoints_hrs = []
for item in timepoints:
    timepoints_datetime = datetime.strptime(item, "%m/%d/%Y %I:%M:%S %p")
    time_delta = timepoints_datetime - initial_time
    timepoints_hrs.append(time_delta.total_seconds() / 3600)

data_dict_processed = process_data_to_mean(data_dict)
#print(data_dict_processed)

#writing data
writing_dict = {}
for key, value in data_dict_processed.items():
    writing_dict[key] = [value[0].tolist(),value[1].tolist(),value[2].tolist()]

with open("output_diya.csv", "w", newline="") as output:
    writer = csv.writer(output, delimiter=",")
    writer.writerow(["timepoints",timepoints_hrs])
    for key, value in writing_dict.items():
        writer.writerow([key,json.dumps(value)])

#plotting
#data_dict_processed.pop("media_met+") #removing unwanted lines
#data_dict_processed.pop("media_met-")

plot_selection = {
    "cells":["JBL137","JBL36"],
    "color":["green","red"],
    "intensity":[2.8,2.3,1.8,1.3,0.8,0.3,0]
}

plot_exclude = {
    "cells":[],
    "color":[],
    "intensity":[]
}

#od600
fig, axs = plt.subplots()
for key, value in data_dict_processed.items():
    key_name_met = key.split("_")
    if key_name_met[0] in plot_exclude["cells"] or key_name_met[1] in plot_exclude["color"] or key_name_met[2] in plot_exclude["intensity"]:
        continue
    else:
        axs.errorbar(timepoints_hrs, value[0][0],
                    yerr = value[1][0], capsize = 2.0,
                    label = key,
                    color = color_map[key_name_met[1]],
                    marker = marker_map[key_name_met[0]], markersize = 3.0,
                    linestyle = linestyle_map[key_name_met[0]], linewidth = 1.0,
                    alpha = alpha_map[key_name_met[2]])

axs.set_xlabel("Time (hrs)")
axs.set_ylabel("OD600")
axs.set_title(f"Cell OD600")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#gfp
fig, axs = plt.subplots()
for key, value in data_dict_processed.items():
    key_name_met = key.split("_")
    if key_name_met[0] in plot_exclude["cells"] or key_name_met[1] in plot_exclude["color"] or key_name_met[2] in plot_exclude["intensity"]:
        continue
    else:
        axs.errorbar(timepoints_hrs, value[0][1],
                    yerr = value[1][1], capsize = 2.0,
                    label = key,
                    color = color_map[key_name_met[1]],
                    marker = marker_map[key_name_met[0]], markersize = 3.0,
                    linestyle = linestyle_map[key_name_met[0]], linewidth = 1.0,
                    alpha = alpha_map[key_name_met[2]])

axs.set_xlabel("Time (hrs)")
axs.set_ylabel("Fluorescence (au)")
axs.set_title(f"Cell GFP")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#gfp/od600
fig, axs = plt.subplots()
for key, value in data_dict_processed.items():
    key_name_met = key.split("_")
    if key_name_met[0] in plot_exclude["cells"] or key_name_met[1] in plot_exclude["color"] or key_name_met[2] in plot_exclude["intensity"]:
        continue
    else:
        axs.errorbar(timepoints_hrs, value[0][2],
                    yerr = value[1][2], capsize = 2.0,
                    label = key,
                    color = color_map[key_name_met[1]],
                    marker = marker_map[key_name_met[0]], markersize = 3.0,
                    linestyle = linestyle_map[key_name_met[0]], linewidth = 1.0,
                    alpha = alpha_map[key_name_met[2]])

#axs.set_xlim(left=4)
#axs.set_ylim(top=2000, bottom=750)
axs.set_xlabel("Time (hrs)")
axs.set_ylabel("Fluorescence (au)/OD600")
axs.set_title(f"Cell GFP/OD600")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()