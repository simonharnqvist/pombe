#!/usr/bin/env Rscript

if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse') # tidying data, ggplot, etc
if (!require('pls')) install.packages('pls'); library('pls') # PLSR modelling
if (!require('Metrics')) install.packages('Metrics'); library('Metrics') # Model performance metrics
if (!require('magrittr')) install.packages('magrittr'); library('magrittr') # Pipes

# Load imputed data
load("../data/temp_data/imputed.Rda")

# Train test split - must ensure that each gene is either in train or test, not both
train_size <- floor(0.75 * nrow(imputed))
index <- sample(c(1:nrow(imputed)), size = train_size, replace = FALSE)
train <- imputed[index, ]
test <- imputed[-index, ]

# Remove redundant first column
train <- select(train, -1)
test <- select(test, -1)

# Remove systematic ID
train <- select(train, -Systematic_ID)
test <- select(test, -Systematic_ID)

# Train PCR model
mod <- plsr(mean.phylop ~ ., data = train, scale = TRUE, center = TRUE, validation = "CV")

# Get optimum number of components (lowest RMSE); this appears to be 8 components
summary(mod) # this can't be saved, must be processed interactively
  
# Retrain with optimum number of components
mod8 <- plsr(mean.phylop ~ ., data = train, scale = TRUE, center = TRUE, validation = "CV", ncomp = 8)

# Predict on test data and get performance metrics
predictions <- predict(mod8, test)
rmse(test$mean.phylop, predictions) # 0.1652 - almost no overfitting

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

# Get variance explained by component
cumul_variance_expl_by_comp <- c(28.933, 36.010, 38.265, 39.34, 39.76, 40.04, 40.34, 40.53) # from summary(mod)
variance_by_component <- variance_expl_per_comp(cumul_variance_expl_by_comp)
write.csv(variance_by_component, "../data/final_data/variance_by_component.csv") # save variance per component to csv

# Get absolute loadings by component
absolute_loadings <- abs_loadings_per_comp(mod8)
absolute_loadings %>% add_column(var = colnames(train)) %>% write.csv(., "../data/final_data/absolute_loadings_per_variable.csv")

# Get percent loadings by component
percent_loadings <- perc_loadings_per_comp(absolute_loadings)
percent_loadings %>% add_column(var = colnames(train)) %>% write.csv(., "../data/final_data/percent_loadings_per_variable.csv")

# Get total variance explained by variable
variance_by_variable <- perc_var_expl_by_variable(mod8, percent_loadings, variance_by_component)
variance_by_variable %>% add_column(var = colnames(train)) %>% write.csv(., "../data/final_data/variance_explained_per_variable.csv")