"""
Hierarchical cluster analysis on conditions perturbed with continual
Ca(OH)2 perturbations with multiple rates, experiment 13
"""

import os
from pathlib import Path

from processing_scripts_formose import (
    comp_info,
    config_file,
    data_report,
    plotting_functions,
)

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage


# Get the path to the data
config = config_file.load_config("./info_files/dir_data.csv")

data_folder = Path(config["dir_extendend_data"])
output_folder = Path(config["output_dir"])
os.makedirs(output_folder / "cluster_analysis", exist_ok=True)

# A list of experiment codes which will be used in the analysis.
experiments = [
    "EXP013",
]

list_exp, exp_condition = config_file.load_exp_info(
    "./info_files/list_exp.csv", experiments
)

# Load compound info
c_info = comp_info.information("./info_files")

# clustering parameters
metric = "correlation"
algorithm = "average"

for exp in experiments:
    # Set up for loading files
    working_path = data_folder / exp / "Analysed_data"
    file_name = working_path / f"{exp}_Data.csv"

    data = data_report.data_report(file=file_name)
    compounds = [comp.split("/")[0] for comp in data.data]

    # Hierarchical clustering analysis
    # Create the distance matrix
    data_as_array = data.to_numpy()
    dist_matrix = pdist(data_as_array, metric)

    # Create the linkage matrix
    Z = linkage(dist_matrix, algorithm, metric, optimal_ordering=False)

    # Plot the dendrogram
    dendrogram = plotting_functions.dendrogram_plot(
        Z,
        [c_info.ind[c_info.SMILES.index(comp)] for comp in compounds],
        f"{str(output_folder)}/cluster_analysis/{exp}_dendrogram",
    )
