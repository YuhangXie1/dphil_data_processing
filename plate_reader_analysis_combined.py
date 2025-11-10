import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import copy
import json
from datetime import datetime

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
            if row_num == 25:
                timepoints_timeformat = row[1]
            if row_num in location_of_data_in_sheet_od600[1]:
                od600_map.loc[row[0]] = np.array(row[1::], dtype=float)
            if row_num in location_of_data_in_sheet_gfp[1]:
                gfp_map.loc[row[0]] = np.array(row[1::], dtype=float)

    return od600_map, gfp_map, timepoints_timeformat

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
filepaths = ["output_plate_dark.csv","output_plate_green.csv","output_plate_red.csv"]
mapping = ["dark","green","red"]
data_dict = {} #[green_JBL131_met+ : [[mean],[std],[CV]]] where [mean] = [OD600, gfp, OD600/gfp]
timepoints = {}
for file_index, file in enumerate(filepaths):
    data_csv = csv.reader(open(file, "r"))
    for rownum, row in enumerate(data_csv):
        if row[0] == "timepoints":
            timepoints[mapping[file_index]+"_"+"timepoints"] = json.loads(row[1])
        else:
            data_dict[mapping[file_index]+"_"+row[0]] = json.loads(row[1])


#plot by strain
plots_by_strain = {} #label (JBL137) : [original key (dark_JBL131_met+), value] where value is the same as data_dict value
for key, value in data_dict.items():
    key_name_strip = key.split("_") #[dark, JBL131, met+]

    if key_name_strip[1] not in plots_by_strain.keys():
        plots_by_strain[key_name_strip[1]] = [[key, value]]
    else:
        plots_by_strain[key_name_strip[1]].append([key, value])

for key, value in plots_by_strain.items():
    key_name_split = value[0].split("_") #[dark, JBL131, met+]
    fig, axs = plt.subplots()
    
    for entry in value[1]:
        axs.errorbar(timepoints[], value[0][0],
                    yerr = value[1][0], capsize = 2.0,
                    label = key,
                    color = color_map[key_name_met[0]],
                    marker = marker_map[key_name_met[1]], markersize = 3.0,
                    linestyle = linestyle_map[key_name_met[1]], linewidth = 1.0)

axs.set_xlabel("Time (hrs)")
axs.set_ylabel("OD600")
axs.set_title(f"Cell OD600 in Dark")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#gfp
fig, axs = plt.subplots()
for key, value in data_dict_processed.items():
    key_name_met = key.split("_")
    axs.errorbar(timepoints_hrs, value[0][1],
                 yerr = value[1][1], capsize = 2.0,
                 label = key,
                 color = color_map[key_name_met[0]],
                 marker = marker_map[key_name_met[1]], markersize = 3.0,
                 linestyle = linestyle_map[key_name_met[1]], linewidth = 1.0)

axs.set_xlabel("Time (hrs)")
axs.set_ylabel("Fluorescence (au)")
axs.set_title(f"Cell GFP in Dark")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()

#gfp/od600
fig, axs = plt.subplots()
for key, value in data_dict_processed.items():
    key_name_met = key.split("_")
    axs.errorbar(timepoints_hrs, value[0][2],
                 yerr = value[1][2], capsize = 2.0,
                 label = key,
                 color = color_map[key_name_met[0]],
                 marker = marker_map[key_name_met[1]], markersize = 3.0,
                 linestyle = linestyle_map[key_name_met[1]], linewidth = 1.0)

axs.set_xlabel("Time (hrs)")
axs.set_ylabel("Fluorescence (au)/OD600")
axs.set_title(f"Cell GFP/OD600 in Dark")
axs.legend(bbox_to_anchor=(1.0, 1.05))

fig.tight_layout()
plt.show()