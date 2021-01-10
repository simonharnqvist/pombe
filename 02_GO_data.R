if (!require(tidyverse)) install.packages('tidyverse'); library(tidyverse) # data manipulation, plotting, etc
if (!require(mltools)) install.packages('mltools'); library(mltools) # one hot encoding
if (!require(data.table)) install.packages('data.table'); library(data.table) # one hot encoding


### GET GENE LIST
# Load Grech data
grech <- read.csv("../data/manually_retrieved_data/hermesgenedata.20180511.txt.gz", sep = "\t")

# Get only gene list
genes <- grech$gene

# Save gene list as tsv
write.table(genes, "../data/genelist.tsv", quote = FALSE, row.names = FALSE, col.names = FALSE)


### READ GO ANNOTATION DATA

# Function to read annotation files
read_annotation_data <- function(file) {
  GO <- read.table(file, sep = "\t") %>%
    select(V2, V9) %>% rename(Term = V2, Genes = V9) %>% .[-1, ]
}

process_GO <- read_annotation_data("../data/GO_slim_process_22112020.txt")
function_GO <- read_annotation_data("../data/GO_slim_function_22112020.txt")
component_GO <- read_annotation_data("../data/GO_slim_component_22112020.txt")

# Split rows to get GO/gene combinations
process_GO_split <- separate_rows(process_GO, Genes) %>% mutate_all(as.factor)
function_GO_split <- separate_rows(function_GO, Genes) %>% mutate_all(as.factor)
component_GO_split <- separate_rows(component_GO, Genes) %>% mutate_all(as.factor)

rm(process_GO, function_GO, component_GO)



### ONE HOT ENCODDE GO ANNOTATION DATA

# Function to one hot encode per gene
one_hot_encode <- function(df, GO_type) {
  
  one_hot_GO <- one_hot(as.data.table(df), "Term") %>%
    as.data.frame %>% group_by(Genes) %>% summarise_all(sum)
  
  colnames(one_hot_GO) <- gsub("Term", GO_type, colnames(one_hot_GO))
  
  return(one_hot_GO)
}

process_onehot <- one_hot_encode(process_GO_split, "Process")
function_onehot <- one_hot_encode(function_GO_split, "Function")
component_onehot <- one_hot_encode(component_GO_split, "Component")

# ENCODE CHROMOSOME (sometimes manually is easier!)
chr <- grech %>% select(gene, chr)
chr$chrI[chr$chr == "I"] <- 1
chr$chrII[chr$chr == "II"] <- 1
chr$chrIII[chr$chr == "III"] <- 1
chr <- chr %>% select(-chr)

# MERGE
GO_onehot <- merge(process_onehot, function_onehot, by = "Genes", all = TRUE) %>% merge(., component_onehot, by = "Genes", all = TRUE) %>%
  merge(., chr, by.x = "Genes", by.y = "gene")

# Set all NAs to 0 (i.e. "not")
GO_onehot[is.na(GO_onehot)] <- 0

# Remove 'none' and then columns with only zeros
GO_onehot <- GO_onehot %>% filter(!Genes == "none") %>% .[, colSums(. != 0) > 0]

# Save file
save(GO_onehot, file = "../data/GO_one_hot_encoded.Rda")
