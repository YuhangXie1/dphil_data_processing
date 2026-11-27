"""
plate_reader_extraction.py
==========================

Step 1 of the plate-reader pipeline.

Reads one Excel export per timepoint (e.g. ``t0.xlsx``, ``t3.xlsx`` ... ``t24.xlsx``)
from a 96-well plate reader and reshapes it into tidy, well-by-time CSV tables
that ``plate_reader_analysis_*.py`` can load.

Pipeline overview
-----------------
    t0.xlsx ... t24.xlsx                      (raw plate reader exports, one per read)
        |
        v   plate_reader_extraction.py        (this script)
        |
    <folder>_<extra>_extracted_<label>.csv    (one table per measurement label,
                                               e.g. OD600, GFP 395nm, GFP 488nm)
    <folder>_<extra>_extracted_combined.csv   (all labels side by side, 2-row header)
        |
        v   plate_reader_analysis_*.py        (+ experiment_input.xlsx plate map)
        |
    figures/*.png

Output table layout
-------------------
Per-label CSV (e.g. ``..._extracted_OD600.csv``)::

            , 19/11/2026 11:47:31, 19/11/2026 15:24:55, ...
        A1  , 0.1286             , 0.1338             , ...
        ...
        H12 , ...

    * rows    = the 96 wells, A1 ... H12
    * columns = the "Start Time" of each read, taken from the Excel file

Combined CSV: the per-label tables concatenated horizontally with a two-level
column header (``measurement``, ``timestamp``) and the index named ``well``.

Supported instruments
---------------------
``"downstairs"``: label rows look like ``"Label: OD600"`` in column A and the
    start time is in column B of a row whose column A contains ``"Start Time:"``.
``"upstairs"``:   label rows have ``"Name"`` in column A and the label in
    column B; the start time is in column E of the ``"Start Time"`` row.

In both formats, plate data rows start with the row letter (``A``-``H``) in
column A, followed by the 12 well values in columns B-M.

Usage
-----
1. Edit the USER INPUT block below (instrument, folder, filename tag, file list).
2. Run from the directory that *contains* ``data_folder``::

       python plate_reader_extraction.py

Requires: numpy, pandas, python-calamine (Excel engine).
"""

import numpy as np
import matplotlib.pyplot as plt  # TIDY: unused in this script, can be removed
import pandas as pd
from os import path


### USER INPUT ###

# Which plate reader produced the Excel files: "downstairs" or "upstairs".
# The two instruments lay out their exports differently (see module docstring).
chosen_plate_reader = "downstairs"

# Folder containing the raw Excel files. Output CSVs are written here too.
data_folder = "26-11-19 Kirill M9 exp 2"

# Extra tag inserted into every output filename, e.g. the analyst's name.
save_filename_extra = "diya"

# Excel files to parse, in chronological order. Each file = one plate read.
# Comment out any file you want to skip (e.g. a failed or duplicate read).
data_files = ["t0.xlsx",
            "t3.xlsx",
            "t6.xlsx",
            "t8.xlsx",
            "t12.xlsx",
            "t15.xlsx",
            "t24.xlsx",
            #"t28.xlsx",
            ]

#data_files = ["incubator t24.xlsx",]

##############


def main():
    """
    Run the full extraction: parse every file in ``data_files``, then write
    one CSV per measurement label plus one combined CSV into ``data_folder``.

    Output filenames::

        <data_folder>/<data_folder>_<save_filename_extra>_extracted_<label>.csv
        <data_folder>/<data_folder>_<save_filename_extra>_extracted_combined.csv
    """
    # Output filenames start with the folder (experiment) name
    save_filename_starter = data_folder

    # --- 1. Extract data from every timepoint file ---------------------------
    # Each file adds one new timestamp column to every label's table in the
    # global dictionary ``all_plate_table_df``.
    for filepath in data_files:
        actual_filepath = path.join(data_folder, filepath)
        if chosen_plate_reader == "downstairs":
            extracted_data = read_data_from_excel_downstairs(actual_filepath)
        elif chosen_plate_reader == "upstairs":
            extracted_data = read_data_from_excel_upstairs(actual_filepath)
        # TIDY: no `else` branch - a typo in chosen_plate_reader gives a
        # confusing NameError on the next line. Consider raising ValueError.
        place_data_into_table(extracted_data)

    # --- 2. Save one CSV per measurement label --------------------------------
    data_tables = {}
    for key, table in all_plate_table_df.items():
        # Drop the placeholder column "0" that the template table was created
        # with. Note: this returns a new table; all_plate_table_df is unchanged.
        new_table = table.drop(columns = "0")
        new_table.to_csv(f"{data_folder}/{save_filename_starter}_{save_filename_extra}_extracted_{key}.csv")
        data_tables[key] = new_table

    # --- 3. Save all labels side by side in one combined CSV -----------------
    combined_table = combine_data(data_tables)
    combined_table.to_csv(f"{data_folder}/{save_filename_starter}_{save_filename_extra}_extracted_combined.csv")


### FUNCTIONS ###

def read_data_from_excel_downstairs(filepath):
    """
    Parse one Excel export from the *downstairs* plate reader.

    Scans the first sheet row by row, keeping track of the most recent
    measurement label (``"Label: <name>"``) and start time (``"Start Time:"``).
    Every row whose first cell is a plate row letter (A-H) is stored under the
    current label.

    On the first call only, also creates an empty 96-well table in the global
    ``all_plate_table_df`` for every label found.

    :param filepath: path to the ``.xlsx`` file for one timepoint
    :type filepath: str
    :return: dict mapping label -> list of ``[start_time, row_letter, values]``,
        where ``values`` holds the 12 readings for that plate row. Example::

            {"OD600": [[t, "A", array([...12 values...])],
                       [t, "B", array([...])], ...],
             "GFP 395nm": [...]}
    :rtype: dict[str, list]
    """
    # header=None: read every row as data, so we can search for keywords ourselves
    data = pd.read_excel(filepath, 0, header=None, engine="calamine")
    all_data = {}
    for row_num, row in enumerate(data.values):
        # "Label: OD600" -> "OD600". Remains in effect until the next label row.
        if "Label:" in str(row[0]):
            data_name = row[0].split(":")[-1].strip()

        # Start time of the current measurement block (used as the column header)
        if "Start Time:" in str(row[0]):
            data_time = row[1]

        # A plate data row: first cell is "A".."H", next 12 cells are wells 1-12
        if str(row[0]) in key_rows:
            if data_name not in all_data.keys():
                all_data[data_name] = [[data_time, row[0], row[1:13]]]
            else:
                all_data[data_name].append([data_time, row[0], row[1:13]])

    # First file only: create an empty 96-well table for each label found.
    # TIDY: labels that appear only in later files will not get a table
    # (KeyError in place_data_into_table).
    if not all_plate_table_df:
        for key in all_data.keys():
            all_plate_table_df[key] = plate_table_df.copy(deep=True)

    return all_data


def read_data_from_excel_upstairs(filepath):
    """
    Parse one Excel export from the *upstairs* plate reader.

    Identical logic to :func:`read_data_from_excel_downstairs`, but the label is
    found in column B of a row containing ``"Name"`` and the start time in
    column E of the ``"Start Time"`` row.

    :param filepath: path to the ``.xlsx`` file for one timepoint
    :type filepath: str
    :return: dict mapping label -> list of ``[start_time, row_letter, values]``
    :rtype: dict[str, list]
    """
    data = pd.read_excel(filepath, 0, header=None, engine="calamine")
    all_data = {}
    for row_num, row in enumerate(data.values):
        # Measurement label is in the cell to the right of "Name"
        if "Name" in str(row[0]):
            data_name = row[1]
            #print(data_name)

        # Start time is in column E on this instrument
        if "Start Time" in str(row[0]):
            data_time = row[4]

        # A plate data row: first cell is "A".."H", next 12 cells are wells 1-12
        if str(row[0]) in key_rows:
            if data_name not in all_data.keys():
                all_data[data_name] = [[data_time, row[0], row[1:13]]]
            else:
                all_data[data_name].append([data_time, row[0], row[1:13]])

    # First file only: create an empty 96-well table for each label found
    if not all_plate_table_df:
        for key in all_data.keys():
            all_plate_table_df[key] = plate_table_df.copy(deep=True)

    return all_data


def place_data_into_table(extracted_data):
    """
    Add one timepoint's readings to the global per-label tables.

    For each label (e.g. ``"OD600"``), flattens the eight row entries into a
    ``{well: value}`` dictionary (``"A1"``, ``"A2"`` ... ``"H12"``) and writes
    it into ``all_plate_table_df[label]`` as a new column whose header is the
    start time of that read.

    Tables are modified in place: rows are wells, columns are timestamps.

    :param extracted_data: output of one ``read_data_from_excel_*`` call
    :type extracted_data: dict[str, list]
    :return: the (mutated) global ``all_plate_table_df``
    :rtype: dict[str, pandas.DataFrame]
    """
    for name, data_entry in extracted_data.items():
        relevant_table = all_plate_table_df[name]  # table for this label, e.g. OD600

        # Build {well_id: value}, e.g. {"A1": 0.128, "A2": 0.124, ...}
        temp_dict = {}
        for row in data_entry:
            row_letter = row[1]
            for index in range(0,len(row[2])):
                temp_dict[row_letter + str(index+1)] = row[2][index]

        # New column header = start time of this label's first row (data_entry[0][0]).
        # .map() aligns values to the well index (A1..H12).
        relevant_table[data_entry[0][0]] = relevant_table.index.map(temp_dict)

    return all_plate_table_df


def combine_data(data_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Concatenate the per-label tables side by side into one wide table.

    The dictionary keys become the top column level, giving a two-level column
    header ``(measurement, timestamp)``, e.g. ``("OD600", "19/11/2026 11:47:31")``.
    This is the format ``load_data`` in the analysis script expects.

    :param data_tables: dict mapping label -> well x timestamp table
    :return: combined table, index named ``"well"``
    """
    combined_df = pd.concat(
        data_tables,
        axis=1,
        names=["measurement", "timestamp"],
    )

    combined_df.index.name = "well"
    return combined_df


### PLATE MAP INITIALISATION (module-level globals) ###

# Plate row letters and column numbers for a standard 96-well plate
key_rows = ["A", "B", "C", "D", "E", "F", "G", "H"]
key_columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

# All 96 well IDs in reading order: A1, A2, ..., A12, B1, ..., H12
key_wells = [str(row)+str(col) for row in key_rows for col in key_columns]

plate_map = {key : np.nan for key in key_wells}  # TIDY: unused, can be removed

# Template table: 96 wells x one placeholder column "0" (dropped before saving)
plate_table_df = pd.DataFrame(0.0, index = key_wells, columns = ["0"])

# Filled in as files are read: {label: well x timestamp DataFrame}
all_plate_table_df = {}

if __name__ == "__main__":
    main()
