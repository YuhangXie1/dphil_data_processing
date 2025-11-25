import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import copy
import json
from datetime import datetime
from os import listdir, path
import re
import xarray as xr

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

#plate map initialisation
key_rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
key_columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
key_wells = [str(row)+str(col) for row in key_rows for col in key_columns]


plate_map = {key : np.nan for key in key_wells}
#print(plate_map)

filepath = "25-11-19_diya_4/25-11-18_diya2_dose_curve_low_gain_t0.xlsx"





#iterates through excel file and extracts all plate data under their respective labels
data = pd.read_excel(filepath, 0, header=None, engine="calamine")
all_data = {}
for row_num, row in enumerate(data.values):
    if "Label" in str(row[0]):
        data_name = row[0].removeprefix("Label:") #finds label
    if "Start Time:" in str(row[0]):
        data_time = row[1] #finds time stamp
    if str(row[0]) in key_rows:
        if data_name not in all_data.keys():
            all_data[data_name] = [data_time, row[0], row[1:13]] #extracting all data into that label
        else:
            all_data[data_name].append([data_time, row[0], row[1:13]])
    
#print(all_data)

#generates a plate map dict for each label, and then fills it with the correct well data
plate_maps = {}
for label, data_value in all_data.items():
    temp_plate_map = {}
    for row_letter in key_rows:
        if row_letter in data_value[1]:
            for index in range(0,len(data_value[2])):
                temp_plate_map[row_letter+str(index+1)] = data_value[2][index]

    print(temp_plate_map)    
    plate_maps[label] = temp_plate_map










#gathering data from each time point and putting into dict
data_dict = {} #name : [od600, gfp]
timepoints = []

    
od600_map, gfp_map, timepoints_timeformat = load_plate_into_df(data)
timepoints.append(timepoints_timeformat)



#working out the time points
initial_time = datetime.strptime(timepoints[0], "%m/%d/%Y %I:%M:%S %p")
timepoints_hrs = []
for item in timepoints:
    timepoints_datetime = datetime.strptime(item, "%m/%d/%Y %I:%M:%S %p")
    time_delta = timepoints_datetime - initial_time
    timepoints_hrs.append(time_delta.total_seconds() / 3600)

data_dict_processed = process_data_to_mean(data_dict)
#print(data_dict_processed)




"""
data = xr.DataArray(np.zeros(2,3),
                    dims = ("condition",
                            "repeat",
                            "channel",
                            "color",
                            "intensity",
                            "time"),

                    coords = {"condition":[],
                              "repeat":[1,2,3],
                              "channel":["OD600","GFP-395","GFP-475"],
                              "color":["green","red","dark"],
                              "intensity":[],
                              "time":[],
                              })
print(data)
print(data.sel(x="a"))
"""
"""
#color and linestyles
color_map = {"green": "green",
             "red": "red",
}

linestyle_map = {"JBL137": "solid",
                 "JBL36": "solid",
                 "JBL001":"dashed",
}

alpha_map = {"2.8": 1,
                 "2.0": 0.8,
                 "1.5": 0.7,
                 "1.0": 0.5,
                 "0.5": 0.4,
                 "0.3": 0.2,
}

marker_map = {"JBL137": "o",
              "JBL36": "^",
              "JBL001":"o",
}

#filepaths = ["25-10-07_growth/25-10-07 Dark t0.xlsx","25-10-07_growth/25-10-07 Dark t2.xlsx","25-10-07_growth/25-10-07 Dark t4.xlsx","25-10-07_growth/25-10-07 Dark t19.xlsx","25-10-07_growth/25-10-07 Dark t22.xlsx"]
#filepaths = ["25-10-07_growth/25-10-07 Green t0.xlsx","25-10-07_growth/25-10-07 Green t2.xlsx","25-10-07_growth/25-10-07 Green t4.xlsx","25-10-07_growth/25-10-07 Green t19.xlsx","25-10-07_growth/25-10-07 Green t22.xlsx"]
#filepaths = ["25-10-28_diya/25-10-28_diya_t0.xlsx","25-10-28_diya/25-10-28_diya_t1.xlsx","25-10-28_diya/25-10-28_diya_t1.5.xlsx","25-10-28_diya/25-10-28_diya_t2.xlsx","25-10-28_diya/25-10-28_diya_t3.xlsx"]
filepath = "25-10-30_diya_2"
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
    "cells":["JBL137","JBL001"],
    "color":[],
    "intensity":[]
}

#od600
fig, axs = plt.subplots()
for key, value in data_dict_processed.items():
    key_name_met = key.split("_")
    print(key_name_met)
    if key == "media":
        axs.errorbar(timepoints_hrs, value[0][0],
                    yerr = value[1][0], capsize = 2.0,
                    label = key,
                    color = "grey")
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
    if key == "media":
        axs.errorbar(timepoints_hrs, value[0][1],
                    yerr = value[1][1], capsize = 2.0,
                    label = key,
                    color = "grey")
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
    if key == "media":
        axs.errorbar(timepoints_hrs, value[0][2],
                    yerr = value[1][2], capsize = 2.0,
                    label = key,
                    color = "grey")
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
"""