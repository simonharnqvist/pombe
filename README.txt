01_retrieve_data.R - gets data from PomBase/STRING via FTP download; also provides a link to Grech data but this is broken (found downloaded data in ../data); only use this script for use with most recent data. For reproducibility, please use data downloaded in ../data
02_GO_data.R - formats gene ontology data
03_format_data.R - formats data to correct format for further applications
04_network_parameters.R - performs network analysis and calculate network centrality measures for each protein (gene)
05_impute_missing_data.R - imputes all missing data with the missForest algorithm
06_variable_importance_estimation.R - trains a principal components regression model; estimates influence of each variable (i.e. variance explained per variable)
07_statistical_relationships.R - outputs a variety of statistical data; correlations between each continuous variable and evolution rate; statistical comparisons between categories
08_plots.R - as the name might suggest, generates plots