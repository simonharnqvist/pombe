library(mltools)
library(tidyverse)
library(doParallel)
library(missForest)
library(Metrics)
library(reticulate)
library(pls)
library(plsVarSel)
library(Metrics)

# Function to read annotation files
read_annotation_data <- function(file) {
  GO <- read.table(file, sep = "\t") %>%
    select(V2, V9) %>% dplyr::rename(Term = V2, Genes = V9) %>% .[-1, ]
}


# Function to one hot encode per gene
one_hot_encode <- function(df, GO_type) {
  
  one_hot_GO <- mltools::one_hot(data.table::as.data.table(df), "Term") %>%
    as.data.frame %>% dplyr::group_by(Genes) %>% dplyr::summarise_all(sum)
  
  colnames(one_hot_GO) <- gsub("Term", GO_type, colnames(one_hot_GO))
  
  return(one_hot_GO)
}

# Function to create graph of interactome
create_interactome <- function(df, conf_threshold) {
  g <- df %>% filter(combined_score >= conf_threshold) %>% igraph::graph_from_data_frame()
  return(g)
}

# Function to get network centralities from graph g
get_centralities <- function(g) {
  nodes <- data.frame(gene = (igraph::V(g)$name)) %>%
    add_column(deg_centr = igraph::degree(g), betw_centr = igraph::betweenness(g),
               closeness_centr = igraph::closeness(g),
               eigen_centr = igraph::eigen_centrality(g)$vector) # add centralities
  return(nodes)
}


# Missing forest imputation
impute_missing <- function(df, n_cores, random_seed, cols_to_ignore) {
  
  impute_df <- df[, !colnames(df) %in% cols_to_ignore] # drop cols
  doParallel::registerDoParallel(cores=n_cores) # for parallelisation
  set.seed(random_seed)
  
  # Impute with missForest
  mf <- missForest::missForest(impute_df, parallelize = "variables")
  
  # Put the removed columns back
  imputed <- cbind(mf$ximp, df[, colnames(df) %in% cols_to_ignore])
  
  
}

# Train test split
train_test_split <- function(df, test_size) {
  
  spec <- c(train = (1-test_size), test = test_size)
  
  # Splitter object
  splitter <- sample(cut(
    seq(nrow(df)),
    nrow(df)*cumsum(c(0,spec)),
    labels = names(spec)
  ))
  
  res <- split(df, splitter)
  
  return(res)
  
  
}

evaluate_mvr <- function(mod, test, y_col, name) {
  # Get optimum number of components
  opt_comps <- pls::selectNcomp(mod)
  
  # CHECK PERFORMANCE ON TEST SET
  y_pred <- predict(mod, test, ncomp = opt_comps)
  rmse_opt <- Metrics::rmse(test[[y_col]], y_pred)
  
  # Pseudo R-squared
  var_expl <- pls::R2(mod, estimate = "test", newdata = test, ncomp = opt_comps)
  
  mod_perf <- data.frame(model = character(), y = character(), ncomps = character(),
                         RMSE = numeric(), var_expl = numeric()) %>%
    dplyr::add_row(model = name, y = y_col, ncomps = as.character(opt_comps), RMSE = rmse_opt,
                   var_expl = var_expl$val[2])
  
  return(mod_perf)
  
}

evaluate_rf <- function(mod, test, y_col, name = "RF") {
  X_cols_test <- test %>% select(-y_col)
  y_pred_rf <- predict(mod, X_cols_test)
  
  # Metrics
  test_y <- test[[y_col]]
  rf_var_exp <- 1 - (mse(test_y, y_pred_rf) / var(test_y))
  rf_rmse <- rmse(test_y, y_pred_rf)
  
  mod_perf <- data.frame(model = character(), y = character(), ncomps = character(),
                         RMSE = numeric(), var_expl = numeric()) %>%
    dplyr::add_row(model = name, y = y_col, ncomps = "None", RMSE = rf_rmse,
                   var_expl = rf_var_exp)
  
  return(mod_perf)
  
}

get_VIP_scores <- function(pls_mod, comps) {
  var_imp <- plsVarSel::VIP(mod, comps) %>% data.frame(t(.)) %>%
    select(1) %>% rownames_to_column()
  colnames(var_imp) <- c("var", "VIP_score")
  return(var_imp)
}

get_percentage_variance_explained <- function(pcr_mod, comps)



fit_random_forest <- function(train, test, y_col) {
  X_cols <- train %>% select(-y_col)
  rf <- randomForest::randomForest(data = train, x = X_cols, y = y_col,
                     ntree = 500, mtry = floor(ncol(X_cols)/3))
  
  X_cols_test <- test %>% select(-y_col)
  y_pred_rf <- predict(rf, X_cols_test)
  
  # Metrics
  rf_var_exp <- 1 - mse(test_y, y_pred_rf) / var(test_y)
  rf_rmse <- rmse(test_y, y_pred_rf)
  
  mod_perf <- data.frame(model = character(), ncomps = numeric(),
                         RMSE = numeric(), var_expl = numeric()) %>%
    dplyr::add_row(model = "RF", ncomps = "None", RMSE = rf_rmse,
            var_expl = rf_var_exp)
  
}

create_variable_categories <- function(df, listOfColNames, varGroupName) {
  category <- df[df$var %in% listOfColNames, ] %>% 
    dplyr::add_column(var_group = varGroupName)
  return(category)
}


get_variance_per_component <- function(mod, opt_comps, cumul_var_by_comp) {
  
  var_by_comp = c()
  cumul_var_by_comp <- as.numeric(cumul_var_by_comp)
  for (i in 2:opt_comps ){
    var_by_comp[i] <- as.numeric(cumul_var_by_comp[i]) - as.numeric(cumul_var_by_comp[i-1])
    var_by_comp[1] <- cumul_var_by_comp[1]
  }
  print(var_by_comp)
}

get_projections <- function(mod, opt_comps) {
  p = mod$projection
  
  p <- p[,1:opt_comps]
  
  for (i in 1:opt_comps){
    p[,i] <- p[,i]^2/t(p[,i]) * p[,i]
  }
  print(p)
}

