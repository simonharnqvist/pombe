if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse') # tidying data, ggplot, etc
if (!require('igraph')) install.packages('igraph'); library('igraph') # network analysis

# Read in Rda interactions data and 'combined' df
load("../data/combined.Rda")
load("../data/interactors.Rda")

# Function to create graph
create_interactome <- function(df, conf_threshold) {
  g <- df %>% filter(Score >= conf_threshold) %>% graph_from_data_frame()
  return(g)
}

# Function to get network centralities from graph g
get_centralities <- function(g) {
  nodes <- data.frame(gene = (V(g)$name)) %>%
    add_column(deg_centr = degree(g), betw_centr = betweenness(g), closeness_centr = closeness(g), eigen_centr = eigen_centrality(g)$vector) # add centralities
  return(nodes)
}

# Apply functions with confidence threshold 400 (moderate confidence)
interactome_data <- interactors %>% create_interactome(., 400) %>% get_centralities()

# Reformat to correct gene names; remove species identifier "4896" and unnecessary trailing ".1" (last two characters)
interactome_data$gene <- interactome_data$gene %>% str_remove_all(., "4896.") %>% str_sub(., end = -3)

# Merge with combined dataframe; save to RDA
combined_full <- merge(combined, interactome_data, by.x = "Systematic_ID", by.y = "gene", all = TRUE)
save(combined_full, file = "../data/combined_full.Rda")


