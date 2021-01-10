if (!require(tidyverse)) install.packages('tidyverse'); library(tidyverse) # data manipulation, plotting, etc

# Read in dataframes
peptide_stats <- read.csv("../data/peptide_stats.tsv", sep = "\t")
grech <- read.csv("../data/hermesgenedata.20180511.txt.gz", sep = "\t") %>% 
  rename(., Systematic_ID = gene) %>% select(-starts_with("hmm"))
interactors <- read.csv("../data/interactors.txt.gz", sep = "", skip = 1, header = FALSE) %>%
  rename(., Protein1 = V1, Protein2 = V2, Score = V3)
aa_composition <- read.csv("../data/aa_composition.tsv", sep = "\t")
load("../data/GO_one_hot_encoded.Rda")

# Fix incorrect formatting in Grech; filter out all rows where all except gene and essential are missing
columns_NA_in_grech <- grech %>% select(-Systematic_ID, -essential) %>% colnames(.)
grech_fixed <- grech[rowSums(is.na(grech[,columns_NA_in_grech])) != length(columns_NA_in_grech),] %>%
  select(-chr) # remove column no longer needed

# Convert AA composition to fraction of total length of protein
aa_matrix <- aa_composition %>% column_to_rownames(var = "Systematic_ID") %>% as.matrix %>% t
aa_scaled <- scale(aa_matrix, center = FALSE, scale = colSums(aa_matrix)) %>% t %>% data.frame %>% rownames_to_column(., var = "Systematic_ID")

# Perform merge each in turn (save all genes from Grech but none of the others)
combined <- merge(grech_fixed, peptide_stats, by = "Systematic_ID", all.x = TRUE) %>%
  merge(., aa_scaled, by = "Systematic_ID", all.x = TRUE) %>%
  merge(., GO_onehot, by.x = "Systematic_ID", by.y = "Genes", all.x = TRUE) %>%
  select(-sum.mRNA.cpc, -start, -end) # Remove column that duplicates RPKM information; start/end could seriously overfit 


save(combined, file = "../data/combined.Rda")

# Save interactions data for network analysis
save(interactors, file = "../data/interactors.Rda")

