#!/usr/bin/env Rscript

library(tidyverse)
library(jsonlite)
library(seqinr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Usage: hbv_resistance.R <path_to_json_folder>", call.=FALSE)
}

# TODO:
# Determine the genotype (is in the json_files) and then select the correct reference for blast
# Needs to be done in the blast process

rtDomain_fasta <- args[1] # rtDomain_fasta <- "/home/jonr/Prosjekter/learning_nextflow/Data/HBV_references/DPOL_HBVD1_RT_domain.fasta"
json_dir       <- "json_files/" # json_dir <- "/home/jonr/delete/"
json_files     <- list.files(path = json_dir,
                             pattern = "json",
                             full.names = TRUE)

# First create a vector of the reference sequence with each amino acid as an element
rtDomain_seq <- unlist(
  # split the sequence between every character and convert to a matrix
  str_split(
    read.fasta(rtDomain_fasta, seqtype = "AA", as.string = TRUE)$RT_domain[1],
    pattern = "",
    simplify = TRUE
  )
)

# Then read pairwise alignment from blastx
# Loop through the json files
for (i in 1:length(json_files)) {
  blast_json <- jsonlite::fromJSON(json_files[i])
  # Get the hsps and convert to tibble
  query_aligned_seq <- unlist(
    str_split(
      # Extract the query hsp
      blast_json[["BlastOutput2"]][["report"]][["results"]][["bl2seq"]][[1]][["hits"]][[1]][["hsps"]][[1]][["qseq"]],
      pattern = "",
      simplify = TRUE
    )
  )
  
  # Find the differences between the aligned sequence and the reference
  tmp <- as_tibble(rtDomain_seq) %>% 
    pivot_longer(everything()) %>% 
    mutate("position" = str_remove(name, "V")) %>% 
    rename("ref_aa" = "value") %>% 
    select(position, ref_aa)
  
  tmp2 <- as_tibble(query_aligned_seq) %>% 
    pivot_longer(everything()) %>% 
    mutate("position" = str_remove(name, "V")) %>% 
    rename("sample_aa" = "value") %>% 
    select(position, sample_aa) %>% 
    # Convert everything to upper characters
    mutate(sample_aa = toupper(sample_aa))
  
  # join the two and create mutations
  df <- left_join(tmp, tmp2) %>% 
    # Create a column indicating when the reference and the sample is not identical
    mutate(mut = case_when(
      ref_aa != sample_aa ~ "mutated"
    )) %>% 
    filter(mut == "mutated") %>% 
    # Create the mutation names
    unite("Mutation", c(ref_aa, position, sample_aa), remove = FALSE, sep = "") %>% 
    select(-mut)
  
  # Write to file
  # Extract the sample name from the json file name
  sampleName <- str_split(basename(json_files[i]), pattern = "\\.", simplify = TRUE)[1,1]
  outfile <- paste0(sampleName, "_HBV_rt_domain_resistance.tsv")
  write_tsv(df, outfile)
}

      
