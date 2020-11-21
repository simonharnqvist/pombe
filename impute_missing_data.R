#!/usr/bin/env Rscript

if (!require('missForest')) install.packages('missForest'); library('missForest') # missing data imputation
if (!require('doParallel')) install.packages('doParallel'); library('doParallel') # parallel processing
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse') # data manipulation, plotting, etc
if (!require('magrittr')) install.packages('magrittr'); library('magrittr') # pipes

# Load dataset
load("../data/processed_data/combined_full.Rda")

# Change datatype of factors to factor
should_be_factor <- c("essential", "GO", "chr")
combined_full[,should_be_factor] %<>% lapply(function(col) as.factor(col))

# Remove gene name and GO from imputation; GO would simply make it take too long, gene name prevents the algorithm from working
without_GO <- combined_full %>% select(-GO) %>% unique()
impute_df <- without_GO %>% select(-Systematic_ID)

# Run imputation (with parallel processing)
registerDoParallel(cores = 3)
set.seed(42)
mf <- missForest(impute_df, parallelize = "variables")

# Combine gene names back with imputed data, then add GO information back on
imputed_without_GO <- cbind(without_GO$Systematic_ID, mf$ximp) %>% rename(without_GO$Systematic_ID)
ID_and_GO <- select(combined_full, Systematic_ID, GO) 
imputed <- list(ID_and_GO, imputed_without_GO) %>% reduce(left_join, by = "Systematic_ID")

# Save imputed for further analysis
save(imputed, file = "../data/processed_data/imputed.Rda")

