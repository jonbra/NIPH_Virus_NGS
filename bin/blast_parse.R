#!/usr/bin/env Rscript
  
library(tidyverse)
library(seqinr)

# Skip this in nextflow
setwd("/home/jonr/Prosjekter/learning_nextflow/HCV_results/5_blast/")

# Read all blast files into a list
tmp <- list.files(full.names = FALSE, pattern = "out$") %>%
  set_names() %>% # To get the file names as the names of the list elements
  map(read_tsv, col_names = FALSE)

# Reduce the list to a single dataframe. Keep the filenames in column 1
df <- bind_rows(tmp, .id = "sampleName")

# Add column for the genotype
df <- df %>% 
  # Fix the sampleName
  mutate(sampleName = str_remove(sampleName, "_blast.out")) %>% 
  separate(X2, into = c("genotype", NA), remove = FALSE) %>% 
  # Count the occurence of each blast subject (i.e. genotype)
  # To know which of the subtypes are most common
  group_by(sampleName) %>% 
  add_count(X2) %>% 
  ungroup()

# Reduce to info on genotypes per group
geno_list <- df %>% 
  group_by(sampleName) %>% 
  distinct(genotype) %>% 
  group_split()

# Write genotypes to file
for (i in 1:length(geno_list)) {
  write_tsv(geno_list[[i]], file = paste0(geno_list[[i]][["sampleName"]][1], "_genotypes.tsv"))
}

# What is the most common subtype per genotype?
subtypes <- df %>% 
  group_by(sampleName, genotype) %>% 
  slice(max(n, n = 1)) %>% 
  distinct(X2)
#subtypes <- scaffolds %>% 
#  group_by(genotype) %>% 
#  slice_max(n, n = 1) %>% 
#  distinct(X2)

# Read the reference fasta file
fasta <- read.fasta(file = "../../Data/Blast_db/HCVgenosubtypes_8.5.19_clean.fa")

# Write out the name of each subtype and the corresponding fasta file
for (i in 1:nrow(subtypes)) {
  write_tsv(subtypes[i, 1], file = paste0(subtypes[i, 1], subtypes[i, 2], ".txt"), col_names = FALSE)
  write.fasta(sequences = fasta[[subtypes$X2[i]]], names = subtypes$X2[i], file.out = paste0(subtypes[i, 1], "_", subtypes[i, 2], "_ref.fa"))
} 

# Divide the scaffolds into the different genotypes

# NB! Dette er spades scaffolds som blir lest?? Blir ikke riktig å splitte slik?
# Split scaffolds per subtype
# Read scaffolds fasta
scaffolds_fa <- read.fasta(file = "${scaffolds_fasta}")

# Split scaffolds per genotype
split <- scaffolds %>% 
  group_split(genotype)

# Write one fasta per genotype
for (i in 1:length(split)){
  # Store the scaffold names and corresponding genotype
  tmp <- split[[i]] %>% select(X1, genotype)
  # Subset the scaffolds
  geno_fa <- scaffolds_fa[tmp\$X1]
  # Write to file
  write.fasta(sequences = geno_fa, names = names(geno_fa), file.out = paste0("${sampleName}_", tmp\$genotype[1], "_scaffolds.fa"))
}

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")