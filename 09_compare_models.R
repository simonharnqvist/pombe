if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse) # tidying data, ggplot, etc
if (!require('randomForest')) install.packages('randomForest'); library(randomForest)
if (!require('pls')) install.packages('pls'); library(pls)
if (!require('Metrics')) install.packages('Metrics'); library(Metrics)

# Get data
load("../data/train.Rda")
load("../data/test.Rda")
mod_perf <- read_csv("../data/model_performance.csv")

# RANDOM FOREST
X <- select(train, -mean.phylop)
y <- unlist(select(train, mean.phylop))
rf <- randomForest(data = train, x = X, y = y, ntree = 500, mtry = floor(ncol(X)/3))

test_X <- select(test, -mean.phylop)
test_y <- unlist(select(test, mean.phylop)) 

y_pred_rf <- predict(rf, test_X)

# Pseudo R-squared
rf_var_exp <- 1 - mse(test_y, y_pred_rf) / var(test_y)

# Get RMSE
rf_rmse <- rmse(test_y, y_pred_rf)

# Add performance to df
mod_perf <- mod_perf %>% add_row(model = "RF", ncomps = NA, RMSE = rf_rmse, var_expl = rf_var_exp)


# PCA REGRESSION
pcr <- pcr(mean.phylop ~ ., data = train, scale = TRUE, center = TRUE, validation = "CV")

# Select optimum components
opt_comps <- selectNcomp(pcr)

# CHECK PERFORMANCE ON TEST SET
y_pred <- predict(pcr, test, ncomp = opt_comps)
rmse_opt <- rmse(test$mean.phylop, y_pred)

# R2
var_expl <- R2(pcr, estimate = "test", newdata = test, ncomp = opt_comps)

# Export performance data
mod_perf <- mod_perf %>%
  add_row(model = "PCR", ncomps = opt_comps, RMSE = rmse_opt, var_expl = var_expl$val[2]) %>%
  write_csv(., "../data/model_performance.csv")

