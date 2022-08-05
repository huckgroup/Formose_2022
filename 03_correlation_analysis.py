"""
 Correlation analysis on difference between means of different sizes (corresponding to time interval)
     e.g. for the 150 s time interval: the difference between means of the samples from t: 0 - 150 s
     and t: 150 - 300 s, t: 30 - 180 s and t: 180 - 330 s, etc., or for the 60 s time interval:the
     difference between means of the samples from t: 0 - 60 s and t: 60 - 120 s, t: 60 - 120 s and
     120 - 180 s, etc.

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
os.makedirs(output_folder / "correlation_analysis", exist_ok=True)


# A list of experiment codes which will be used in the analysis.
experiments = ["EXP013"]

list_exp, exp_condition = config_file.load_exp_info(
    "./info_files/list_exp.csv", experiments
)

# Load compound info
c_info = comp_info.information("./info_files")

l = list(zip(c_info.ind, c_info.SMILES))
time_intervals = [150, 120, 90, 60, 30]  # in seconds
sample_time = 30  # in seconds


for exp in exp_condition:
    working_path = data_folder / exp / "Analysed_data"
    file_name = working_path / f"{exp}_Data.csv"

    data = data_report.data_report(file=file_name)
    compounds = [*data.data]

    flow = data.conditions["NaOH_flow/ µl/h"]  # each step is 1 second
    flow_time = data.conditions["flow_profile_time/ s"]

    # create dictionarry with flow values at data points of sample times ###
    flow_values = dict()
    flow_values["data_points"] = []
    for x in data.series_values:
        flow_values["data_points"].append(flow[int(x)])

    # different differentials on mean bins for time intervals in 'time_intervals' ###
    d_data = data_analysis_functions.differential_means(
        data.data, time_intervals, sample_time, l
    )
    d_flow = data_analysis_functions.differential_means(
        flow_values, time_intervals, sample_time, [("no_ind", "data_points")]
    )

    # Pearson correlation analysis for different time scales ###
    corr = data_analysis_functions.correlation(d_data, d_flow, time_intervals)

    indexes = []
    for a, b in l:
        for x in [*data.data]:
            if b in x:
                indexes.append(a)

    # file writer   ['correlation_value','hex_colour']
    file_writers.write_corr_csv(
        corr,
        time_intervals,
        indexes,
        filename=output_folder / "correlation_analysis" / f"correlation_analysis.csv",
    )
