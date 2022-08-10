# Analysis Software Accompanying the Publication "Trajectory Dependent Compositional Outcomes of a Prebiotic Reaction Network"

## Installation Instructions

Using the anaconda package manager to create a virtual environment and install
the following packages from the command line:

```
conda env create fnoise-env
conda activate fnoise-env
conda install scipy seaborn
pip install statannotations
```

Change the paths `info_files/dir_data.csv`. `dir_extendend_data` is the path to
the directory containing the extended data on your computer. `output_dir` is
the path to the directory in which you would like analysis output files (plots
and .csv results) stored on your computer.

These scripts can be run from the command line from their parent directory:

```
python <script_name>.py
```

## Description of scripts

### 01_composition_analysis.py

Run this script to generate violin plots for series of data in Figure 2 and
Figures S6-S11.

### 02_time_trace_hierarchical_clustering.py

This script performs hierarchical clustering of traces within each experiment.
The output is a dendrogram plot for EXP013, Figure 3b.

### 03_correlation_analysis.py

This script performs the time-interval correlation analysis (see Materials and
Methods). It outputs the results of the correlation analysis as a .csv file,
which is used to create Figure 3c.

### 04_compositional_shift.py

This script generates plots indicating how groups of compounds respond
collectively to applied perturbations, used to create Figure 5a.

