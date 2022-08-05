"""
Composition shift from steady state on data sets with combined experiments:
    2 - 3 with respect to 1,
    5 - 6 with respect to 4,
    8 - 9 with respect to 7,
    and 11 - 12 with respect to 10,
    normalized between -1 and 1 for each compound
"""

import os
from pathlib import Path

from processing_scripts_formose import (
    comp_info,
    data_analysis_functions,
    config_file,
    data_report,
    file_writers,
)

# Get the path to the data
config = config_file.load_config("./info_files/dir_data.csv")

data_folder = Path(config["dir_extendend_data"])
output_folder = Path(config["output_dir"])
os.makedirs(output_folder / "compositional_shift", exist_ok=True)

# A list of experiment codes which will be used in the analysis.
experiment_sets = [
    ["EXP004", "EXP005", "EXP006"],  # 20 mM
    ["EXP001", "EXP002", "EXP003"],  # 50 mM
    ["EXP007", "EXP008", "EXP009"],  # 100 mM
    ["EXP010", "EXP012", "EXP011"],  # 50 mM freq variation
]

experiment_list = [
    "EXP005",
    "EXP006",
    "EXP002",
    "EXP003",
    "EXP008",
    "EXP009",
    "EXP012",
    "EXP011",
]


# Load compound info
c_info = comp_info.information("./info_files")

compound_colours = comp_info.load_colours_dict("./info_files/compound_information.csv")

dic_diff = dict()
for x in c_info.ind:
    dic_diff[x] = []

compound_numbering = list(zip(c_info.ind, c_info.SMILES))

for c, set in enumerate(experiment_sets, 0):
    current_set = []  # store for the data in each series
    for exp in set:

        # Load data
        working_path = data_folder / exp / "Analysed_data"

        file_name = working_path / f"{exp}_Data.csv"

        data = data_report.data_report(file=file_name)

        current_set.append(data)

    # calculate relative difference perturbed state from steady state
    dic_diff = data_analysis_functions.difference_average(
        current_set[0], current_set[1], compound_numbering, dic_diff
    )
    dic_diff = data_analysis_functions.difference_average(
        current_set[0], current_set[2], compound_numbering, dic_diff
    )

# normalize relative difference from steady state for each compound, between -1 and 1, inf values are replaced with 1
dic_rel_diff = data_analysis_functions.normalized_difference(
    dic_diff, compound_numbering
)

file_writers.write_rel_diff_csv(
    dic_rel_diff,
    experiment_list,
    filename=output_folder
    / "compositional_shift"
    / f"relative_concentration_differences.csv",
)
