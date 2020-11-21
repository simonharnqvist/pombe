#!/usr/bin/env Rscript

if (!require(RCurl)) {install.packages("RCurl") library(RCurl)}
if (!require(tidyverse)) {install.packages("tidyverse") library(tidyverse)}

# Set wd
dir.create("../pombase_data")
setwd("../pombase_data")

# Note: conservation and general genomics data from Grech et al 2019 is retrieved manually from 
#download.file("https://ndownloader.figshare.com/files/11450096", destfile = "grech.txt.gz") # - this is causing issues when downloaded automatically (file naming issue?)

# Get protein data from PomBase
peptide_stats <- "ftp://ftp.pombase.org/pombe/Protein_data/PeptideStats.tsv"
download.file(peptide_stats, destfile = "peptide_stats.tsv")

aa_composition <- "ftp://ftp.pombase.org/pombe/Protein_data/aa_composition.tsv"
download.file(aa_composition, destfile = "aa_composition.tsv")

# Get GO data from PomBase
GO <- "ftp://ftp.pombase.org/pombe/annotations/Gene_ontology/gene_association.pombase.gz"
download.file(GO, destfile = "gene_ontology.gz")

# Get physical interactors from STRING
interactors <- "https://stringdb-static.org/download/protein.links.v11.0/4896.protein.links.v11.0.txt.gz"
download.file(interactors, destfile = "interactors.txt.gz")
