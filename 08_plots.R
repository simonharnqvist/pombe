if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse) # data manipulation, plotting, etc
if (!require('Rmisc')) install.packages('Rmisc'); library(Rmisc) # summary stats for plots
if (!require('viridis')) install.packages('viridis'); library(viridis) # colours for plots

# Load dataframes
load("../data/VIP.Rda")
load("../data/correlations.Rda")
load("../data/continuous_wilcoxon.Rda")

# PLOT: VIP BY VARIABLE GROUP

# Get statistics for error bars  
summary_varimps <- summarySE(variable_importance, measurevar = "VIP_score", groupvars = "var_group")

# Plot
var_imp_plot <- variable_importance %>%
  ggplot(aes(x = reorder(var_group, VIP_score), y = VIP_score)) +
  geom_point(position = "jitter", size = 1.5, colour = "cyan3") + 
  coord_flip() +
  geom_errorbar(data = summary_varimps, aes(x = var_group, ymin = VIP_score, ymax = VIP_score, width = 0.3), col = "black", size = 1.2) +
  geom_errorbar(data = summary_varimps, aes(x = var_group, ymin = VIP_score - se, ymax = VIP_score + se, width = 0.3), col = "black", size = 1.2) +
  theme_classic() +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15)) +
  xlab("Variable group") + ylab("Variable importance in projection")

# Save plot
png(file = "../plots/vip.png", width = 600, height = 400)
var_imp_plot
dev.off()


# PLOT: CORRELATIONS

varnames <- c("Degree centrality", "Gene expression", "Betweeness centrality",
              "Alanine", "Glycine", "Valine", "Gene length", "Protein length", "Protein mass",
              "Arginine", "Aspartic acid", "Glutamic acid", "Protein expr", 
              "Lysine", "Glutamine", "Methionine", "Isoleucine", "Histidine", "Threonine",
              "Tryptophan", "Proline", "pI", "Charge",
              "CAI", "Tyrosine", 
              "Leucine", "Cysteine", "Phenylalanine",
              "Asparagine", "Serine", "Knockout fitness")

correlations$Significance[correlations$p_adj < 0.05] <- "Significant"
correlations$Significance[correlations$p_adj >= 0.05] <- "Not significant"

correlations_plot <- correlations %>% filter(x_var != "mean.phylop") %>%
  ggplot(aes(x = reorder(x_var, cor), y = cor, col = Significance)) +
  scale_x_discrete(labels = rev(varnames)) +
  geom_point(size = 2) + 
  geom_segment(aes(x = reorder(x_var, -cor), xend = reorder(x_var, -cor), y = 0, yend = cor), size = 1) +
  geom_hline(yintercept = 0, size = 1.2, linetype = "dashed") +
  theme_classic() +
  coord_flip() +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15)) +
  xlab("Variable") + ylab("Correlation with constraint")

# Save
png(file = "../plots/correlations.png", width = 600, height = 500)
correlations_plot
dev.off()



# PLOT: MEAN CHANGE IN PHYLOP WHEN PRESENT

summary_diff <- summarySE(mw_df, measurevar = "mean_diff", groupvars = "var_group")

mw_df$Significance[mw_df$p_adj < 0.05] <- "Significant"
mw_df$Significance[mw_df$p_adj >= 0.05] <- "Not significant"

mean_diff_plot <- ggplot(data = mw_df, aes(x = reorder(var_group, mean_diff), y = mean_diff, col = Significance)) +
  geom_point(position = "jitter", size = 1.5) + 
  geom_errorbar(data = summary_diff, aes(x = var_group, ymin = mean_diff, ymax = mean_diff, width = 0.3), colour = "black", size = 1.2) +
  geom_errorbar(data = summary_diff, aes(x = var_group, ymin = mean_diff - se, ymax = mean_diff + se, width = 0.3), colour = "black", size = 1.2) +
  theme_classic() +
  geom_hline(yintercept = 0, size = 1.2, linetype = "dashed") +
  theme(axis.title = element_text(size = 20), axis.text = element_text(size = 15)) +
  xlab("Variable group") + ylab("Difference in constraint (phyloP)")

# Save
# Save
png(file = "../plots/mean_diff.png", width = 700, height = 400)
mean_diff_plot
dev.off()
