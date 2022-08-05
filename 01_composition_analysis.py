"""
Composition analysis on data sets with combined experiments 1 - 3, 4 - 6, 7 - 9
and 10 - 12
"""

import os
from pathlib import Path

from processing_scripts_formose import (
    comp_info,
    data_analysis_functions,
    config_file,
    data_report,
    file_writers,
    plotting_functions,
)

# Get the path to the data
config = config_file.load_config("./info_files/dir_data.csv")

data_folder = Path(config["dir_extendend_data"])
output_folder = Path(config["output_dir"])
os.makedirs(output_folder / "statistics", exist_ok=True)

exp_name = ["01_20_mM_amp", "02_50_mM_amp", "03_100_mM_amp", "04_50_mM_freq"]
for x in exp_name:
    os.makedirs(output_folder / "violin_plots" / f"{x}", exist_ok=True)

# A list of experiment codes which will be used in the analysis.
experiment_sets = [
    ["EXP004", "EXP005", "EXP006"],  # 20 mM
    ["EXP001", "EXP002", "EXP003"],  # 50 mM
    ["EXP007", "EXP008", "EXP009"],  # 100 mM
    ["EXP010", "EXP012", "EXP011"],  # 50 mM freq variation
]


# These variables can be placed in a separate configuration file.
independent_variable_units = ["σ/ mM", "σ/ mM", "σ/ mM", "rate/ s"]
independent_variables = [
    ["0", "2.89", "5.75"],
    ["0", "2.89", "5.76"],
    ["0", "2.89", "5.77"],
    ["0", "45", "120"],
]


# Load compound info
c_info = comp_info.information("./info_files")
compound_colours = comp_info.load_colours_dict("./info_files/compound_information.csv")
names = dict(zip(c_info.SMILES, c_info.name))
index = dict(zip(c_info.SMILES, c_info.ind))
compound_numbers = list(zip(c_info.ind, c_info.SMILES))

# The indices of each sequence to compare
pair_indices = [(0, 1), (1, 2), (0, 2)]

exp_idx = 0
for c, set in enumerate(experiment_sets, 0):
    current_set = []  # store for the data in each series
    exp_n = exp_name[exp_idx]
    for exp in set:

        # Load data
        working_path = data_folder / exp / "Analysed_data"

        file_name = working_path / f"{exp}_Data.csv"

        data = data_report.data_report(file=file_name)

        current_set.append(data)

        # calculate averages and standard deviations
        averages = data_analysis_functions.data_averages(data)
        standard_deviations = data_analysis_functions.data_standard_deviations(data)

        # write the averages and standard deviations to files
        file_writers.write_average_stdev_csv(
            averages,
            standard_deviations,
            index,
            filename=output_folder / "statistics" / f"{exp}_statistics.csv",
        )

    # Calculate p_values
    p_values_1_2 = data_analysis_functions.data_p_values(
        current_set[0], current_set[1], compound_numbers
    )
    p_values_2_3 = data_analysis_functions.data_p_values(
        current_set[1], current_set[2], compound_numbers
    )
    p_values_1_3 = data_analysis_functions.data_p_values(
        current_set[0], current_set[2], compound_numbers
    )

    # Generate violin plots
    plotting_functions.create_series_violin_plots(
        current_set,
        compound_colours=compound_colours,
        series_values=independent_variables[c],
        x_label=independent_variable_units[c],
        filename=str(
            output_folder / "violin_plots" / f"{exp_n}" / f"{exp_n}_violin_plots"
        ),
        pairs=pair_indices,
        p_values=[p_values_1_2, p_values_2_3, p_values_1_3],
        names=names,
        index=index,
    )
    exp_idx += 1
