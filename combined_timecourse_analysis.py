import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import json

#functions
def load_data_assign_plate(filepath):
    data = csv.reader(open(filepath, 'r'))
    #assigning data to plate
    dark_data_dict = {}
    time_points = []

    for row_num, row in enumerate(data):
        if row_num == 3:
            time_points = np.array(row[1::], dtype=np.float32)
        if row[0] in plate_map.keys():
            dark_data_dict[plate_map[row[0]]] = np.array(row[1::], dtype=np.float32)

    return dark_data_dict, time_points

def mean_std_cv(data):
    #data should be array containing 1 or more arrays
    mean = np.mean(data, axis = 0)
    std_dev = np.std(data, axis = 0)
    cv = (std_dev/mean)*100
    return [mean, std_dev, cv]

def process_data_to_mean(filepath):
    dark_data_dict, time_points = load_data_assign_plate(filepath)
    dark_data_dict_repeats = {}
    for key, value in dark_data_dict.items():
        #parsing name
        namemet_repeat = key.rsplit("_",1)
        if namemet_repeat[0] in dark_data_dict_repeats.keys():
            dark_data_dict_repeats[namemet_repeat[0]].append(value)
        else:
            dark_data_dict_repeats[namemet_repeat[0]] = [value]

    dark_data_dict_processed = {}
    for key, value in dark_data_dict_repeats.items():
        #getting mean, std, cv
        mean_std_cv_array = mean_std_cv(value)

        #calculating upper bound
        upper_bound = mean_std_cv_array[0] + mean_std_cv_array[1]
        lower_bound = mean_std_cv_array[0] - mean_std_cv_array[1]
        mean_std_cv_array.append(upper_bound)
        mean_std_cv_array.append(lower_bound)

        #new dictionary with processed data. value = [mean, std, cv, upper bound, lower bound]
        dark_data_dict_processed[key] = mean_std_cv_array

    return dark_data_dict_processed, time_points

### MAIN ###

#plate map
plate_map = {"B1":"media_met+_all",
             "C1":"media_met+_kan",
             "D1":"media_met+_none",
             "B12":"media_met-_all",
             "C12":"media_met-_kan",
             "D12":"media_met-_none",

             "B3":"JBL137_met+_1",
             "B4":"JBL137_met+_2",
             "B5":"JBL137_met+_3",

             "B8":"JBL137_met-_1",
             "B9":"JBL137_met-_2",
             "B10":"JBL137_met-_3",

             "D3":"JBL131_met+_1",
             "D4":"JBL131_met+_2",
             "D5":"JBL131_met+_3",

             "D8":"JBL131_met-_1",
             "D9":"JBL131_met-_2",
             "D10":"JBL131_met-_3",

             "F3":"JBL066_met+_1",
             "F4":"JBL066_met+_2",
             "F5":"JBL066_met+_3",

             "F8":"JBL066_met-_1",
             "F9":"JBL066_met-_2",
             "F10":"JBL066_met-_3",

             "H3":"JBL001_met+_1",
             "H4":"JBL001_met+_2",
             "H5":"JBL001_met+_3",

             "H8":"JBL001_met-_1",
             "H9":"JBL001_met-_2",
             "H10":"JBL001_met-_3",
             }

#color and linestyles
color_map = {"JBL137": "blue",
             "JBL131": "orange",
             "JBL066": "brown",
             "JBL001": "black",
             "media": "yellow",
}

linestyle_map = {"met+": "solid",
                 "met-": "dashed",
}

marker_map = {"met+": "o",
              "met-": "^",
}

#loading time point data
timepoints_data = csv.reader(open("output.csv", "r"))
timepoints_data_dict = {}
for rownum, row in enumerate(timepoints_data):
    timepoints_data_dict[row[0]] = json.loads(row[1])
timepoints_data_timepoints = timepoints_data_dict.pop("timepoints")

#plotting OD600
filepath = "25-10-07 Repeat wellplate growth trials od600.csv"
dark_data_dict_processed_od600, time_points = process_data_to_mean(filepath)
#calculating OD offset between timepoints and timecourse
offset = np.average(dark_data_dict_processed_od600["media_met+"][0]) - np.average(timepoints_data_dict["media_met+"][0][0])
dark_data_dict_processed_od600.pop("media_met+") #removing unwanted lines
dark_data_dict_processed_od600.pop("media_met-")

fig, axs = plt.subplots()
for key, value in dark_data_dict_processed_od600.items():
    key_name_met = key.split("_")
    axs.plot(time_points, value[0], label = key, color = color_map[key_name_met[0]], linestyle = linestyle_map[key_name_met[1]])
    axs.fill_between(time_points, value[3], value[4], alpha = 0.1, color = color_map[key_name_met[0]])

#timepoint data
for key, value in timepoints_data_dict.items():
    key_name_met = key.split("_")
    axs.errorbar(timepoints_data_timepoints, value[0][0]+offset, yerr = value[1][0], label = key, color = color_map[key_name_met[0]], marker = marker_map[key_name_met[1]], linestyle = '')

axs.set_xlabel("Time (hrs)")
axs.set_ylabel("OD600")
axs.set_title(f"Cell OD600 in Dark")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#plotting GFP
filepath = "25-10-07 Repeat wellplate growth trials GFP.csv"
dark_data_dict_processed_gfp, time_points = process_data_to_mean(filepath)
dark_data_dict_processed_gfp.pop("media_met+") #removing unwanted lines
dark_data_dict_processed_gfp.pop("media_met-")
fig, axs = plt.subplots()
for key, value in dark_data_dict_processed_gfp.items():
    key_name_met = key.split("_")
    axs.plot(time_points, value[0], label = key, color = color_map[key_name_met[0]], linestyle = linestyle_map[key_name_met[1]])
    axs.fill_between(time_points, value[3], value[4], alpha = 0.1, color = color_map[key_name_met[0]])

for key, value in timepoints_data_dict.items():
    key_name_met = key.split("_")
    axs.errorbar(timepoints_data_timepoints, value[0][1], yerr = value[1][1], label = key, color = color_map[key_name_met[0]], marker = marker_map[key_name_met[1]], linestyle = '')


axs.set_xlabel("Time (hrs)")
axs.set_ylabel("GFP")
axs.set_title(f"Cell GFP in Dark")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#GFP/OD600
dark_data_gfp_div_od600 = {}
for key, value in dark_data_dict_processed_gfp.items():
    dark_data_gfp_div_od600[key] = [value[i]/dark_data_dict_processed_od600[key][i] for i in range(0,len(value))]

fig, axs = plt.subplots()
for key, value in dark_data_gfp_div_od600.items():
    key_name_met = key.split("_")
    axs.plot(time_points, value[0], label = key, color = color_map[key_name_met[0]], linestyle = linestyle_map[key_name_met[1]])
    axs.fill_between(time_points, value[3], value[4], alpha = 0.1, color = color_map[key_name_met[0]])

for key, value in timepoints_data_dict.items():
    key_name_met = key.split("_")
    axs.errorbar(timepoints_data_timepoints, value[0][2], yerr = value[1][2], label = key, color = color_map[key_name_met[0]], marker = marker_map[key_name_met[1]], linestyle = '')


axs.set_xlabel("Time (hrs)")
axs.set_ylabel("GFP/OD600")
axs.set_title(f"Cell GFP/OD600 in Dark")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()




