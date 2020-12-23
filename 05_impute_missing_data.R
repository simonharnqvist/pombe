#!/usr/bin/env Rscript

if (!require('missForest')) install.packages('missForest'); library('missForest') # missing data imputation
if (!require('doParallel')) install.packages('doParallel'); library('doParallel') # parallel processing
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse') # data manipulation, plotting, etc
if (!require('magrittr')) install.packages('magrittr'); library('magrittr') # pipes

# Load dataset
load("../data/temp_data/combined_full.Rda")

# Change datatype of factors to factor
should_be_factor <- combined_full %>%
  select(essential, starts_with("Process"), starts_with("Function"), starts_with("Component"), starts_with("chr")) %>% colnames
combined_full[,should_be_factor] %<>% lapply(function(col) as.factor(col))

# Remove gene name as this prevents algorithm from working; remove GO and amino acid information in the interest of
# saving time
impute_df <- combined_full %>% select(-Systematic_ID, -starts_with("Function"), -starts_with("Process"),
                                      -starts_with("Component"), -A, -C, -D, -E, -F, -G, -M, -H, -I, -K, -L,
                                      -N, -P, -Q, -R, -S, -T, -V, -W, -Y)

# Run imputation (with parallel processing)
registerDoParallel(cores = 3)
set.seed(42)
mf <- missForest(impute_df, parallelize = "variables")

# Combine gene names back with imputed data
imputed <- cbind(combined_full$Systematic_ID, mf$ximp) %>% cbind(., combined_full %>% select(Systematic_ID, starts_with("Function"), starts_with("Process"), 
                                                                                             starts_with("Component"), A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y))
# Where GO data are missing, impute 0 (absent) - i.e. any remaining NAs = 0
imputed[is.na(imputed)] <- 0 

# Remove redundant column
imputed <- imputed %>%  select(-`combined_full$Systematic_ID`)

# Save imputed for further analysis
save(imputed, file = "../data/temp_data/imputed.Rda")

# Save as CSV
write.csv(imputed, "../data/final_data/imputed.csv")

