#!/usr/bin/env Rscript

if (!require(tidyverse)) install.packages('tidyverse'); library(tidyverse) # data manipulation, plotting, etc
if (!require(RCurl)) install.packages('RCurl'); library(RCurl) # retrieve data via URLs

# Set wd
dir.create("../data/pombase_data")
setwd("../data/pombase_data")

# Get protein data from PomBase
peptide_stats <- "ftp://ftp.pombase.org/pombe/Protein_data/PeptideStats.tsv"
download.file(peptide_stats, destfile = "peptide_stats.tsv")

aa_composition <- "ftp://ftp.pombase.org/pombe/Protein_data/aa_composition.tsv"
download.file(aa_composition, destfile = "aa_composition.tsv")

# Get physical interactors from STRING (featured in PomBase)
interactors <- "https://stringdb-static.org/download/protein.links.v11.0/4896.protein.links.v11.0.txt.gz"
download.file(interactors, destfile = "interactors.txt.gz")
