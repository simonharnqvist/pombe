require(tidyverse)

# Set wd
setwd("../data")

# Read in dataframes
peptide_stats <- read.csv("peptide_stats.tsv", sep = "\t")
gene_ontology <- read.csv("gene_ontology.gz", sep = "\t", skip = 24, header = FALSE) %>% select(V2, V5) %>% rename(., Systematic_ID = V2, GO = V5)
grech <- read.csv("hermesgenedata.20180511.txt.gz", sep = "\t") %>% rename(., Systematic_ID = gene) %>% select(-starts_with("hmm"))
interactors <- read.csv("interactors.txt.gz", sep = "", skip = 1, header = FALSE) %>% rename(., Protein1 = V1, Protein2 = V2, Score = V3)
aa_composition <- read.csv("aa_composition.tsv", sep = "\t")

# Convert AA composition to fraction of total length of protein
aa_matrix <- aa_composition %>% column_to_rownames(var = "Systematic_ID") %>% as.matrix %>% t
aa_scaled <- scale(aa_matrix, center = FALSE, scale = colSums(aa_matrix)) %>% t %>% data.frame %>% rownames_to_column(., var = "Systematic_ID")

# Perform outer join on all except interactions data; save
combined <- list(peptide_stats, gene_ontology, grech, aa_scaled) %>% reduce(left_join, by = "Systematic_ID")
save(combined, file = "../data/processed_data/combined.Rda")

# Save interactions data for network analysis
save(interactors, file = "../data/processed_data/interactors.Rda")

