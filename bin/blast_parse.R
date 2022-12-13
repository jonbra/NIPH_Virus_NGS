#!/usr/bin/env Rscript
  
library(tidyverse)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: blast_parse.R <sampleName> <blast_out> <scaffolds> <references>", call.=FALSE)
}

sampleName <- args[1]
blast_out  <- args[2]
scaffolds  <- args[3]
references <- args[4]

scaf <- read_tsv(blast_out, col_names = FALSE) %>% 
  # Add a column for the genotype
  separate(X2, into = c("genotype", NA), remove = FALSE) %>% 
  # Count the number of each subject
  add_count(X2) 

# Which genotypes are present?
genotypes <- scaf %>% 
  distinct(genotype)

write_tsv(genotypes, paste0(sampleName, ".genotypes.tsv"))

# What is the most common subtype per genotype?
subtypes <- scaf %>% 
  group_by(genotype) %>% 
  slice_max(n, n = 1) %>% 
  distinct(X2)

# Read the reference fasta file
fasta <- read.fasta(file = references)

# Write out the name of each subtype
for (i in 1:nrow(subtypes)) {
  write_tsv(subtypes[i, 1], file = paste0(sampleName, ".", subtypes$X2[i], ".txt"), col_names = FALSE)
  write.fasta(sequences = fasta[[subtypes$X2[i]]], names = subtypes$X2[i], file.out = paste0(sampleName, ".", subtypes$X2[i], "_ref.fa"))
} 

# NB! Dette er spades scaffolds som blir lest?? Blir ikke riktig å splitte slik?
# Split scaffolds per subtype
# Read scaffolds fasta
scaffolds_fa <- read.fasta(file = scaffolds)

# Split scaffolds per genotype
split <- scaf %>% 
  group_split(genotype)

# Write one fasta per genotype
for (i in 1:length(split)){
  # Store the scaffold names and corresponding genotype
  tmp <- split[[i]] %>% select(X1, genotype)
  # Subset the scaffolds
  geno_fa <- scaffolds_fa[tmp$X1]
  # Write to file
  write.fasta(sequences = geno_fa, names = names(geno_fa), file.out = paste0(sampleName, ".", tmp$genotype[1], "_scaffolds.fa"))
}

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")