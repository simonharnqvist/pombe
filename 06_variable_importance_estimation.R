#!/usr/bin/env Rscript

if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse') # tidying data, ggplot, etc
if (!require('pls')) install.packages('pls'); library('pls') # PLSR modelling
if (!require('Metrics')) install.packages('Metrics'); library('Metrics') # Model performance metrics
if (!require('magrittr')) install.packages('magrittr'); library('magrittr') # Pipes
if (!require('plsVarSel')) install.packages('plsVarSel'); library('plsVarSel') # VIP

# Load imputed data
load("../data/imputed.Rda")

# Get df for modelling
modelling_df <- select(imputed, -Systematic_ID)

# Split datasets into train/test/val
spec <- c(train = .8, test = .2)

splitter <- sample(cut(
  seq(nrow(modelling_df)), 
  nrow(modelling_df)*cumsum(c(0,spec)),
  labels = names(spec)
))

res <- split(modelling_df, splitter)

train <- res$train
test <- res$test

save(train, file = "../data/train.Rda")
save(test, file = "../data/test.Rda")

# Train PLS model
mod <- plsr(mean.phylop ~ ., data = train, scale = TRUE, center = TRUE, validation = "CV")

# Get optimum number of components
opt_comps <- selectNcomp(mod)

# CHECK PERFORMANCE ON TEST SET
y_pred <- predict(mod, test, ncomp = opt_comps)
rmse_opt <- rmse(test$mean.phylop, y_pred)

# Pseudo R-squared
var_expl <- R2(mod, estimate = "test", newdata = test, ncomp = opt_comps)

# Export performance
mod_perf <- data.frame(model = character(), ncomps = numeric(), RMSE = numeric(), var_expl = numeric()) %>%
  add_row(model = "PLS", ncomps = opt_comps, RMSE = rmse_opt, var_expl = var_expl$val[2]) %>%
  write_csv(., "../data/model_performance.csv")

# VARIABLE IMPORANCE IN PROJECTION
var_imp <- VIP(mod, opt_comps) %>% data.frame(t(.)) %>% select(1) %>% rownames_to_column()
colnames(var_imp) <- c("var", "VIP_score")

# Assign broad variable categories
process <- var_imp[grep("Process", var_imp$var),] %>% add_column(var_group = "Process")
funct <- var_imp[grep("Function", var_imp$var),] %>% add_column(var_group = "Function")
component <- var_imp[grep("Component", var_imp$var),] %>% add_column(var_group = "Component")
size <- var_imp %>% filter(var == "genelength" | var == "Mass..kDa." | var == "Residues") %>% add_column(var_group = "Size")
importance <- var_imp %>% filter(var == "essential1" | var == "solid.med.fitness") %>% add_column(var_group = "Functional importance")
expression <- var_imp %>% filter(var == "log.phase.RPKM" | var == "sum.protein.cpc")  %>% add_column(var_group = "Expression")
other <- var_imp %>% filter(var == "pI" | var == "Charge" | var == "CAI") %>% add_column(var_group = "Other")
amino_acids <- var_imp %>% filter(nchar(var) == 1) %>% add_column(var_group = "Amino acid content")
network <- var_imp %>% filter(str_detect(var, "centr"))  %>% add_column(var_group = "Network centrality")
chromosome <- var_imp %>% filter(var == "chrI1" | var == "chrII1" | var == "chrIII1")  %>% add_column(var_group = "Chromosome")

# Combine to single dataframe & save CSV file + R data
variable_importance <- rbind(process, funct, component, size, importance, expression, other, amino_acids, network, chromosome)

write_csv(variable_importance, "../data/variable_importance_in_projection.csv")
save(variable_importance, file = "../data/VIP.Rda")

