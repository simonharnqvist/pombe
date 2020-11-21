#!/usr/bin/env Rscript

if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse') # tidying data, ggplot, etc
if (!require('pls')) install.packages('pls'); library('pls') # PLSR modelling
if (!require('Metrics')) install.packages('Metrics'); library('Metrics') # Model performance metrics

# Load imputed data
load("../data/processed_data/imputed.Rda")

# Train test split - must ensure that each gene is either in train or test, not both
without_GO <- imputed %>% select(-GO) %>% distinct()
train_size <- floor(0.75 * nrow(without_GO))
index <- sample(c(1:nrow(without_GO)), size = train_size, replace = FALSE)
train_without_GO <- without_GO[index, ]
test_without_GO <- without_GO[-index, ]

# Combine with GO in each dataset
train <- imputed %>% select(Systematic_ID, GO) %>% list(train_without_GO, .) %>% reduce(left_join, by = "Systematic_ID")
test <- imputed %>% select(Systematic_ID, GO) %>% list(test_without_GO, .) %>% reduce(left_join, by = "Systematic_ID")

# Train PCR model
mod <- pcr(mean.phylop ~ ., data = train, scale = TRUE, center = TRUE, validation = "CV")

# Get optimum number of components (lowest RMSE)

# Retrain with optimum number of components

# Predict on test data and get performance metrics

# Function to get variance explained by each variable by each component
variance_expl_per_comp <- function(cumul_variance_expl_by_comp) {
  var_by_comp <- c()
  var_by_comp[1] <- cumul_variance_expl_by_comp[1]
  
  # Then get variance explained by each component
  for (i in 2:length(cumul_variance_expl_by_comp)) {
    var_by_comp[i] <- cumul_variance_expl_by_comp[i] - cumul_variance_expl_by_comp[i - 1]
  }
  
  return(var_by_comp)
}

# Function to get variable importance (absolute loadings per component)
abs_loadings_per_comp <- function(model) {
  mat <- as.matrix(abs(model$loadings)) # get matrix of absolute loadings
  
  abs_loadings <- data.frame()
  
  # For each variable, get absolute loadings per comp
  for (i in 1:nrow(mat)) {
    infl_per_comp <- t(mat[i,])
    abs_loadings <- rbind(abs_loadings, infl_per_comp)
  }
  
  return(abs_loadings)
  
}

# Function to get percent loadings
perc_loadings_per_comp <- function(abs_loadings) {
  perc_loadings <- data.frame(apply(abs_loadings[1:ncol(abs_loadings)], 2, FUN = function(x){(x / sum(x)) * 100}))
  
  # Turn index column into proper column, allows grouping of variables later
  perc_loadings <- cbind(Variable = rownames(perc_loadings), perc_loadings)
  rownames(perc_loadings) <- 1:nrow(perc_loadings) # change to numeric row names
  
  return(perc_loadings)
}

# Function to get total variance explained by each variable
perc_var_expl_by_variable <- function(model, perc_loadings, var_by_comp) {
  
  # Multiply each column by the variance explained by that PC
  
  plsr_varimps <- data.frame(t(t(as.matrix(perc_loadings[, 2:ncol(perc_loadings)])) * var_by_comp / 100))
  
  # Get total variance explained per variable
  variance_per_variable <- data.frame(var = rownames(model$loadings), variance_expl = rowSums(plsr_varimps))
  
  return(variance_per_variable)               
}

# Apply functions to model

# Collapse factors to get group level data

# Save dataframes as appropriate
# variance_explained_per_variable_per_component
# variance_explained_per_variable
# variance_explained_per_variable_group

