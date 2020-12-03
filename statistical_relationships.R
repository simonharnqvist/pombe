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

# For all other variables, performing Mann-Whitney with Bonferroni correction is best choice
wilcoxon <- function(df, y) {

  # Select factors that have two levels (0, 1)
  factors <- df %>% select(where(is.factor)) %>% select_if(~ nlevels(.) == 2) %>% colnames
  
  # Create dataframe
  mw_df <- data.frame(x_var = character(), p = numeric(), W = numeric(), mean_0 = numeric(), mean_1 = numeric(), sd_0 = numeric(), sd_1 = numeric())
  
  # Perform Wilcoxon for each x variable
  for (x in factors) {
    mww_mod <- wilcox.test(df[[y]] ~ df[[x]])
  
    # Calculate stats per level
    zero <- df %>% subset(.[[x]] == 0)
    one <- df %>% subset(.[[x]] == 1)
    mean_zero <- mean(zero[[y]])
    mean_one <- mean(one[[y]])
    sd_zero <- sd(zero[[y]])
    sd_one <- sd(zero[[y]])
  
   # Add to dataframe
   mw_df <- mw_df %>% add_row(x_var = as.character(x), p = as.numeric(mww_mod$p.value), W = as.numeric(mww_mod$statistic),
                              mean_0 = as.numeric(mean_zero), mean_1 = as.numeric(mean_one),
                              sd_0 = as.numeric(sd_zero), sd_1 = as.numeric(sd_one))
  }

  
  # Bonferroni correction
  mw_df$p_adj <- mw_df$p * length(factors)
}

categorical <- wilcoxon(imputed, "mean.phylop")




