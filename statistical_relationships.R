if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse) # data manipulation, plotting, etc
if (!require('FSA')) install.packages('FSA'); library(FSA) # Dunn test

# Load data
load("../data/temp_data/imputed.Rda")

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


# Function for Dunn-Kruskal-Wallis test
dunn_kw <- function(df, x, y) {
  x_vctr <- df[[x]]
  y_vctr <- df[[y]]
  
  dunn_mod <- dunnTest(x_vctr, y_vctr)
  return(dunn_mod$res)
}

# Chromosome is the only continuous variable with more than two levels
chromosome_kw <- dunn_kw(imputed, "mean.phylop", "chr")

# Perform Wilcoxon with multiple test correction for all other factors

# Select factors that have two levels (0, 1)
factors <- imputed %>% select(where(is.factor)) %>% select_if(~ nlevels(.) == 2) %>% colnames
  
# Create dataframe
mw_df <- data.frame(x_var = character(), p = numeric(), W = numeric())
  
# Perform Wilcoxon for each x variable (would be better if this could be a function)
for (x in factors) {
  mww_mod <- wilcox.test(imputed[["mean.phylop"]] ~ imputed[[x]])
    
  # Add to dataframe
  mw_df <- mw_df %>% add_row(x_var = x, p = mww_mod$p.value, W = mww_mod$statistic)
}

# Bonferroni correction
mw_df$p_adj <- mw_df$p * length(factors)

# Save files
save(correlations, file = "../data/final_data/correlations.Rda")
save(chromosome_kw, file = "../data/final_data/chromosome.Rda")
save(mw_df, file = "../data/final_data/continuous_wilcoxon.Rda")
