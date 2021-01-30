library(tidyverse)
library(ggpubr)
library(viridis)
library(Rmisc)

# Read data
pcors <- read_csv("../data/partial_correlations.csv")
imp <- read_csv("../data/variable_importance_in_projection.csv")

# Set up palette
palette <- viridis(10, option = "C")

colours <- c("AA composition" = palette[1], "Chromosome" = palette[2],
            "Component" = palette[3], "Expression" = palette[4],
            "Function" = palette[5], "Functional importance" = palette[6],
            "Network centrality" = palette[6], "Other" = palette[7],
            "Process" = palette[8], "Size" = palette[9])


# Merge dataframes
imp <- imp %>% mutate(var = str_replace_all(var, "1|`", "")) # remove unwanted characters preventing merger
plot_df <- merge(imp, pcors) %>% select(var, VIP_score, r, "p-adj", group)

# PLOT: PARTIAL CORRELATIONS
plot_df$logp <- -log10(plot_df$`p-adj`)

labs <- c("A", "Betweenness", "C", "CAI", "Chr I", "Chr III", "Closeness", "Cytosolic", "Mitochondrial",
          "D", "Degree centrality", "Eigencentrality", "G", "Gene length", "Gene expression", 
          "Protein mass", "N", "Biosynthetic process", "Small mol metabol", "R", "Protein length", "S",
          "KO fitness", "Protein expression", "V", "W")

cor_plot <- plot_df %>% filter(plot_df$`p-adj` < 0.05) %>% add_column(labels = labs) %>%
  ggtext(., x = "r", y = "logp", color = "group", label = labs, repel = TRUE, label.rectangle = TRUE) +
  theme_classic() +
  labs(x = "Partial correlation (T) with constraint", y = expression("-log"[10]*"(p.adj.)")) +
  theme(legend.position = "none") +
  scale_colour_manual(values = colours)

# Save
ggsave(plot = cor_plot, "cor_plot.png", path = "../plots", width = 185, height = 159, unit = "mm", dpi = 600)

# PLOT: VIP
vip_plot <- plot_df %>% ggplot(., aes(x = reorder(group, -VIP_score), y = VIP_score, color = group)) +
  geom_point(position = "jitter") +
  scale_colour_manual(values = colours) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), legend.title = element_blank()) +
  labs(x = "Variable group", y = "Variable importance in projection")

# Error bars
summary <- summarySE(plot_df, measurevar = "VIP_score", groupvars = "group")
vip_plot <- vip_plot + geom_errorbar(data = summary, aes(ymin = VIP_score, ymax = VIP_score), width = 0.5, size = 2)    

# Save
ggsave(plot = vip_plot, "vip_plot.png", path = "../plots", width = 185, height = 159, unit = "mm", dpi = 600)
                           