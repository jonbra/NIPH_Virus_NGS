#!/usr/bin/env Rscript
  
library(tidyverse)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Usage: blast_parse.R <sampleName> <blast_out> <scaffolds> <references> <agens>", call.=FALSE)
}

sampleName <- args[1]
blast_out  <- args[2]
scaffolds  <- args[3]
references <- args[4]
agens      <- args[5]

if (agens == "HCV") {

  scaf <- read_tsv(blast_out, col_names = FALSE) %>% 
    # Rename columns
    rename("qseqid" = "X1",
      "sseqid" = "X2",
      "pident" = "X3",
      "length" = "X4",
      "mismatch" = "X5",
      "gapopen" = "X6",
      "qstart" = "X7",
      "qend" = "X8",
      "sstart" = "X9",
      "send" = "X10",
      "evalue" = "X11",
      "bitscore" = "X12") %>% 
    # Add a column for the genotype
    separate(sseqid, into = c("genotype", NA), remove = FALSE) %>% 
    # Count the number of each subject
    add_count(sseqid) 
  
  # Write out the reformatted blast result
  write_csv(scaf %>% select(-n), file = paste0(sampleName, "_blast_out.csv"))

  # Which genotypes are present?
  genotypes <- scaf %>% 
    distinct(genotype)

  write_tsv(genotypes, paste0(sampleName, ".genotypes.tsv"))

  # What is the most common subtype per genotype?
  subtypes <- scaf %>% 
    group_by(genotype) %>% 
    slice_max(n, n = 1) %>% 
    distinct(sseqid)

  # Read the reference fasta file
  fasta <- read.fasta(file = references)

  # Write out the name of each subtype
  for (i in 1:nrow(subtypes)) {
    write_tsv(subtypes[i, 1], file = paste0(sampleName, ".", subtypes$sseqid[i], ".txt"), col_names = FALSE)
    write.fasta(sequences = fasta[[subtypes$sseqid[i]]], names = subtypes$sseqid[i], file.out = paste0(sampleName, ".", subtypes$sseqid[i], "_ref.fa"))
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
    tmp <- split[[i]] %>% select(qseqid, genotype)
    # Subset the scaffolds
    geno_fa <- scaffolds_fa[tmp$qseqid]
    # Write to file
    write.fasta(sequences = geno_fa, names = names(geno_fa), file.out = paste0(sampleName, ".", tmp$genotype[1], "_scaffolds.fa"))
  }  else if (agens == "HBV") {
  blast_res <- read_tsv(blast_out, col_names = FALSE) %>% 
    # Rename columns
    rename("qseqid" = "X1",
           "sseqid" = "X2",
           "pident" = "X3",
           "length" = "X4",
           "mismatch" = "X5",
           "gapopen" = "X6",
           "qstart" = "X7",
           "qend" = "X8",
           "sstart" = "X9",
           "send" = "X10",
           "evalue" = "X11",
           "bitscore" = "X12") %>% 
    # Count the number of each Blast subject hit
    add_count(sseqid) 
  
  # Write out the reformatted blast result
  write_csv(blast_res %>% select(-n), file = paste0(sampleName, "_blast_out.csv"))
  
  # Read the reference fasta file
  fasta <- read.fasta(file = references)
  
  # Get a list of names of the blast hits
  hits <- blast_res %>% 
    distinct(sseqid)
  
  for (i in 1:nrow(hits)) {
    # Write a file with the name of each blast hit reference
    write_tsv(hits[i, 1], file = paste0(sampleName, ".", hits$sseqid[i], ".txt"), col_names = FALSE)
    # Write a fasta file for each blast hit reference
    write.fasta(sequences = fasta[[hits$sseqid[i]]], names = hits$sseqid[i], file.out = paste0(sampleName, ".", hits$sseqid[i], "_ref.fa"))
  } 
}

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")
