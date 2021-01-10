if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse) # data manipulation, plotting, etc

# Load data
load("../data/imputed.Rda")

# Function to get Kendall correlations
kendall_correlations <- function(df, y) {
  
  # Select numeric columns
  num_cols <- df %>% select(where(is.numeric))
  
  # Create dataframe with just col names for correlation data
  correlations <- data.frame(x_var = character(), cor = numeric(), p = numeric(), t = numeric())
  
  # Get correlation stats for each x-y pair
  for (x in colnames(num_cols)) {
    cor_mod <- cor.test(df[[x]], df[[y]], method = "kendall")
    correlations <- correlations %>% add_row(x_var = as.character(x), cor = as.numeric(cor_mod$estimate),
                                             p = as.numeric(cor_mod$p.value), t = as.numeric(cor_mod$statistic))
  }
  
  return(correlations)
}

correlations <- kendall_correlations(imputed, "mean.phylop")

# Bonferroni corrections
correlations$p_adj <- correlations$p * length(correlations$x_var)


# Perform Wilcoxon with multiple test correction for all other factors

# Select factors that have two levels (0, 1)
factors <- imputed %>% select(where(is.factor)) %>% select_if(~ nlevels(.) == 2) %>% colnames
  
# Create dataframe
mw_df <- data.frame(x_var = character(), p = numeric(), W = numeric(), mean_diff = numeric())
  
# Perform Wilcoxon for each x variable (would be better if this could be a function)
for (x in factors) {
  mww_mod <- wilcox.test(imputed[["mean.phylop"]] ~ imputed[[x]])
  mean_diff <- mean(imputed[["mean.phylop"]][imputed[[x]] == 1]) - mean(imputed[["mean.phylop"]][imputed[[x]] == 0])
    
  # Add to dataframe
  mw_df <- mw_df %>% add_row(x_var = x, p = mww_mod$p.value, W = mww_mod$statistic, mean_diff = mean_diff)
}

# Bonferroni correction
mw_df$p_adj <- mw_df$p * length(factors)

# Group by variable group
process <- mw_df[grep("Process", mw_df$x_var),] %>% add_column(var_group = "Process")
funct <- mw_df[grep("Function", mw_df$x_var),] %>% add_column(var_group = "Function")
component <- mw_df[grep("Component", mw_df$x_var),] %>% add_column(var_group = "Component")
essential <- mw_df[grep("essential", mw_df$x_var),] %>% add_column(var_group = "Essentiality")
chromosome <- mw_df %>% filter(x_var == "chrI" | x_var == "chrII" | x_var == "chrIII") %>% add_column(var_group = "Chromosome")

mw_df <- rbind(process, funct, component, essential, chromosome)

# Save files
save(correlations, file = "../data/correlations.Rda")
save(mw_df, file = "../data/continuous_wilcoxon.Rda")
write_csv(correlations, "../data/correlations.csv")
write_csv(mw_df, "../data/binary_differences.csv")
