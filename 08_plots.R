library(Rmisc)
library(viridis)
library(ggpubr)
library(tidyverse)

# Read data
pcors <- read_csv("../data/partial_correlations.csv")
imp <- read_csv("../data/variable_importance_in_projection.csv")

# Set up palette
palette <- viridis(10, option = "C")


#############
### FIGURE 1
#############

colours <- c("AA composition" = palette[1], "Chromosome" = palette[3],
            "Component" = palette[5], "Expression" = palette[7],
            "Function" = palette[9], "Functional importance" = palette[2],
            "Network centrality" = palette[4], "Other" = palette[6],
            "Process" = palette[8], "Size" = palette[10])


# Merge dataframes
imp <- imp %>% mutate(var = str_replace_all(var, "1|`", "")) # remove unwanted characters preventing merger
plot_df <- merge(imp, pcors) %>% select(var, VIP_score, r, "p-adj", group)

# PLOT: PARTIAL CORRELATIONS
plot_df$logp <- -log10(plot_df$`p-adj`)

labs <- c("A", "Betweenness", "C", "CAI", "Chr III", "Closeness", "Cytosolic", "Mitochondrial",
          "D", "Degree centrality", "Eigencentrality", "G", "Gene length", "Gene expression", 
          "Protein mass", "N", "Biosynthetic process", "R", "Protein length", "S",
          "KO fitness", "Protein expression", "V", "W")

cor_plot <- plot_df %>% filter(plot_df$`p-adj` < 0.05) %>% add_column(labels = labs) %>%
  ggtext(., x = "r", y = "logp", color = "group", label = labs, repel = TRUE, label.rectangle = TRUE) +
  theme_classic() +
  labs(x = "Partial correlation (T) with constraint", y = expression("-log"[10]*"(p.adj.)")) +
  theme(legend.position = "none") +
  scale_colour_manual(values = colours)

# PLOT: VIP
vip_plot <- plot_df %>% ggplot(., aes(x = reorder(group, -VIP_score), y = VIP_score, color = group)) +
  geom_point(position = "jitter") +
  scale_colour_manual(values = colours) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  labs(x = "Variable group", y = "Variable importance in projection")

# Error bars
summary <- summarySE(plot_df, measurevar = "VIP_score", groupvars = "group")
vip_plot <- vip_plot + geom_errorbar(data = summary, aes(ymin = VIP_score, ymax = VIP_score), width = 0.5, size = 2)    

# Save
annotate_figure(vip_plot, left = "A")
fig1 <- ggarrange(vip_plot, cor_plot, ncol = 1, nrow = 2, common.legend = FALSE, align = "v")
ggsave(plot = fig1, "fig1.png", path = "../plots", width = 130, height = 200, unit = "mm", dpi = 600)   


#############
### FIGURE 2
#############

# STEP PLOT OF 10 LARGEST VIP

# Correct labels

fi_labs <- c("Cellular nitrogen compound metabolism", 
             "Protein-containing complex", 
             "Translation", 
             "Structural component of ribosome",
             "Ribosome",
             "Intracellular",
             "Biosynthetic process",
             "Cellular",
             "Structural molecule activity",
             "RNA binding")


go_plot <- imp %>% filter(var_group == "Process" | var_group == "Component" | var_group == "Function") %>%
  arrange(desc(VIP_score)) %>% slice(1:10) %>%
  ggplot(., aes(x = reorder(var, VIP_score), y = VIP_score, fill = var_group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = 1, label = fi_labs)) + 
  theme_classic() +
  ylab("Variable importance in projection") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), legend.title = element_blank()) +
  scale_fill_manual(values = colours) +
  coord_flip()

# Save
ggsave(plot = go_plot, "fig2.png", path = "../plots", width = 130, height = 100, unit = "mm", dpi = 600)   


  