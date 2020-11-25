#!/usr/bin/env Rscript

if (!require('missForest')) install.packages('missForest'); library('missForest') # missing data imputation
if (!require('doParallel')) install.packages('doParallel'); library('doParallel') # parallel processing
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse') # data manipulation, plotting, etc
if (!require('magrittr')) install.packages('magrittr'); library('magrittr') # pipes

# Load dataset
load("../data/temp_data/combined_full.Rda")

# Change datatype of factors to factor
should_be_factor <- combined_full %>% select(essential, chr, starts_with("Process"), starts_with("Function"), starts_with("Component")) %>% colnames
combined_full[,should_be_factor] %<>% lapply(function(col) as.factor(col))

# Remove gene name as this prevents algorithm from working
impute_df <- combined_full %>% select(-Systematic_ID)

# Run imputation (with parallel processing)
registerDoParallel(cores = 3)
set.seed(42)
mf <- missForest(impute_df, parallelize = "variables")

# Combine gene names back with imputed data
imputed <- cbind(combined_full$Systematic_ID, mf$ximp)


# Save imputed for further analysis
save(imputed, file = "../data/temp_data/imputed.Rda")

