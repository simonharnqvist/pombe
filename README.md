# Determinants of sequence evolution constraints in the protein coding genome of fission yeast
This repository contains the R scripts for the paper "Determinants of the rate of molecular evolution in fission yeast".

## Retrieve data
All data are found on https://figshare.com/projects/Determinants_of_the_rate_of_molecular_evolution_in_fission_yeast/96245. Some data are also automatically retrieved in the first script. The GO data must be processed manually, unfortunately; instructions for this are given in the manuscript.

## Run scripts
Once all data, including GO Slims, are retrieved, the scripts can be run as they are. The scripts are:
* (**01_retrieve_data.R** - retrieves data from databases, but it is better for reproducibility to not run this script and use the data from Figshare; only run this if the latest data are required)
* **02_GO_data.R** - formats GO data (as well as chromosome data) to one-hot encoded format
* **03_format_data.R** - merges the various data sources into a single table, removes redundant columns, and scales amino acid information to account for protein length
* **04_network_parameters.R** - performs basic network analysis and calculates betweenness and degree centralities for each gene
* **05_impute_missing_data.R** - imputation with missForest algorithm; if multiple cores are available this can be parallelised for speed (currently set to 3 cores)
* **06_variable_importance_estimation.R** - trains a PLS model and calculates Variable Importance in Projection (VIP) for each variable
* **07_statistical_relationships.R** - calculates correlations with constraint for each continous variable; difference in mean constraint for one-hot encoded variables
* **08_plots.R** - draws plots
* **09_model_comparisons.R** - compares PLS model with random forest and PCR models
