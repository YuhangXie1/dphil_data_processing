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
#iterates through excel file and extracts all plate data in rows
def read_data_from_excel(filepath):
    data = pd.read_excel(filepath, 0, header=None, engine="calamine")
    all_data = {}
    for row_num, row in enumerate(data.values):
        if "Label:" in str(row[0]):
            data_name = row[0].split(":")[-1].strip() #finds label
            #print(data_name)
        if "Start Time:" in str(row[0]):
            data_time = row[1] #finds time stamp
        if str(row[0]) in key_rows:
            if data_name not in all_data.keys():
                all_data[data_name] = [[data_time, row[0], row[1:13]]] #extracting all data into that label
            else:
                all_data[data_name].append([data_time, row[0], row[1:13]])

    if not all_plate_table_df:
        for key in all_data.keys():
            all_plate_table_df[key] = plate_table_df.copy(deep=True)
    
    return all_data

def place_data_into_table(extracted_data):
    #pandas dataframe for containing data. Rows are plate coordinates (e.g. A1) Columns are organised time (e.g. 0 hr, 2 hr)
    #initiating tables for each type of data
    #print(extracted_data)
    for name, data_entry in extracted_data.items():
        relevant_table = all_plate_table_df[name] #selects table by name (e.g. OD600)

        temp_dict = {} #makes a dictionary where index (e.g. A1) gets assigned the right data point
        for row in data_entry:
            row_letter = row[1]
            for index in range(0,len(row[2])):
                temp_dict[row_letter + str(index+1)] = row[2][index]   

        relevant_table[data_entry[0][0]] = relevant_table.index.map(temp_dict) #inserts data into the right table

    return all_plate_table_df

### MAIN ###

#plate map initialisation
key_rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
key_columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
key_wells = [str(row)+str(col) for row in key_rows for col in key_columns]
plate_map = {key : np.nan for key in key_wells}

plate_table_df = pd.DataFrame(0.0, index = key_wells, columns = ["0"])
all_plate_table_df = {}
#print(plate_map)

folder = "25-11-27_diya_5"
name = "25-11-27_diya_5"
filepaths = ["25-11-27_diya2_t0.xlsx",
             "25-11-27_diya2_t5.xlsx",
             "25-11-27_diya2_t11.xlsx",
             "25-11-27_diya2_t24.xlsx"]
#filepath = "25-11-19_diya_4/25-11-18_diya2_dose_curve_low_gain_t0.xlsx"

#extracting data
for filepath in filepaths:
    actual_filepath = path.join(folder, filepath)
    extracted_data = read_data_from_excel(actual_filepath)
    place_data_into_table(extracted_data)

#printing data
for key, table in all_plate_table_df.items():
    #deleting first column which is place holder. Note does not update the tables in all_plate_table_df
    new_table = table.drop(columns = "0")
    new_table.to_csv(f"{folder}/{name}_extracted_{key}.csv")