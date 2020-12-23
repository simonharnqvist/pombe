#!/usr/bin/env Rscript

if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse') # tidying data, ggplot, etc
if (!require('pls')) install.packages('pls'); library('pls') # PLSR modelling
if (!require('Metrics')) install.packages('Metrics'); library('Metrics') # Model performance metrics
if (!require('magrittr')) install.packages('magrittr'); library('magrittr') # Pipes
if (!require('plsVarSel')) install.packages('plsVarSel'); library('plsVarSel') # VIP

# Load imputed data
load("../data/temp_data/imputed.Rda")

# Train test split - must ensure that each gene is either in train or test, not both
train_size <- floor(0.75 * nrow(imputed))
index <- sample(c(1:nrow(imputed)), size = train_size, replace = FALSE)
train <- imputed[index, ]
test <- imputed[-index, ]

# Remove systematic ID
train <- select(train, -Systematic_ID)
test <- select(test, -Systematic_ID)

# Train PLS model
mod <- plsr(mean.phylop ~ ., data = train, scale = TRUE, center = TRUE, validation = "CV")

# Get optimum number of components (lowest RMSE); this appears to be 4 components
summary(mod) 

# Predict on test data and get performance metrics
y_pred <- predict(mod, test, ncomp = 5)
rmse(test$mean.phylop, y_pred) # 0.165 - no apparent overfitting

#R sq
1 - mse(test$mean.phylop, y_pred) / var(test$mean.phylop) # 0.352

# Get variable importance in projection
var_imp <- VIP(mod, 5) %>% data.frame(t(.)) %>% select(1) %>% rownames_to_column()
colnames(var_imp) <- c("var", "VIP_score")

# Assign broad variable categories
process <- var_imp[grep("Process", var_imp$var),] %>% add_column(var_group = "Process")
funct <- var_imp[grep("Function", var_imp$var),] %>% add_column(var_group = "Function")
component <- var_imp[grep("Component", var_imp$var),] %>% add_column(var_group = "Component")
size <- var_imp %>% filter(var == "genelength" | var == "Mass..kDa." | var == "Residues") %>% add_column(var_group = "Size")
essential_dispens <- var_imp %>% filter(var == "essential1" | var == "solid.med.fitness") %>% add_column(var_group = "Essentiality/dispensability")
expression <- var_imp %>% filter(var == "log.phase.RPKM" | var == "sum.protein.cpc")  %>% add_column(var_group = "Expression")
other <- var_imp %>% filter(var == "pI" | var == "Charge" | var == "CAI") %>% add_column(var_group = "Other")
amino_acids <- var_imp %>% filter(nchar(var) == 1) %>% add_column(var_group = "Amino acid content")
network <- var_imp %>% filter(str_detect(var, "centr"))  %>% add_column(var_group = "Network centrality")
chromosome <- var_imp %>% filter(var == "chrI1" | var == "chrII1" | var == "chrIII1")  %>% add_column(var_group = "Chromosome")

# Combine to single dataframe & save CSV file + R data
variable_importance <- rbind(process, funct, component, size, essential_dispens, expression, other, amino_acids, network, chromosome)

write_csv(variable_importance, "../data/final_data/variable_importance_in_projection.csv")
save(variable_importance, file = "../data/temp_data/VIP.Rda")

