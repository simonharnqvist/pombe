if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse) # tidying data, ggplot, etc
if (!require('randomForest')) install.packages('randomForest'); library(randomForest)
if (!require('pls')) install.packages('pls'); library(pls)
if (!require('Metrics')) install.packages('Metrics'); library(Metrics)


# Get data
load("../data/temp_data/imputed.Rda")

# TRAIN/TEST SPLIT
train_size <- floor(0.75 * nrow(imputed))
index <- sample(c(1:nrow(imputed)), size = train_size, replace = FALSE)
train <- imputed[index, ]
test <- imputed[-index, ]

train <- select(train, -Systematic_ID)
test <- select(test, -Systematic_ID)

# RANDOM FOREST
X <- select(train, -mean.phylop)
y <- unlist(select(train, mean.phylop))
rf <- randomForest(data = train, x = X, y = y)

test_X <- select(test, -mean.phylop)
test_y <- unlist(select(test, mean.phylop)) 

y_pred_rf <- predict(rf, test_X)

# Get rsq
1 - mse(test_y, y_pred_rf) / var(test_y) # 0.467

# Get RMSE
rmse(test_y, y_pred_rf)


# PCA regression
pcr <- pcr(mean.phylop ~ ., data = train, scale = TRUE, center = TRUE, validation = "CV")

summary(pcr)

# Lowest RMSE
y_pred_1 <- predict(pcr, test, ncomp = 139)
rmse(test$mean.phylop, y_pred_1) # 0.171
1 - mse(test$mean.phylop, y_pred_1) / var(test$mean.phylop) # 0.342

# Save ncomps as PLS
y_pred_2 <- predict(pcr, test, ncomp = 5)
rmse(test$mean.phylop, y_pred_2) # 0.184
1 - mse(test$mean.phylop, y_pred_2) / var(test$mean.phylop) # 0.237

# 90% of available variance explanced (0.9 * max var)
0.9 * 38.73 # 34.9 -> 93 comps
y_pred_3 <- predict(pcr, test, ncomp = 93)
rmse(test$mean.phylop, y_pred_3) # 0.172
1 - mse(test$mean.phylop, y_pred_3) / var(test$mean.phylop) # 0.330



