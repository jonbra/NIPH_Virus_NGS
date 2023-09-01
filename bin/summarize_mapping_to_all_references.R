#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Usage: summarize_mapping_to_all_references.R <idxstats file>", call.=FALSE)
}

idxstats   <- args[1]
sampleName <- args[2]

# Read the summary of the first mapping
df <- read_tsv(idxstats, col_names = FALSE) %>% 
  # Identify genotype and subtype
  separate(X1, into = c("Subtype", "Reference"), sep = "_", remove = FALSE) 

# Count number of reads per subtype
summary <- df %>% 
  group_by(Subtype) %>% 
  summarise(reads = sum(X3)) %>% 
  arrange(desc(reads))

# Find major reference to use
major_tmp <- summary$Subtype[1] 
major_ref <- df %>% 
  filter(Subtype == major_tmp) %>% 
  # Choose the reference with most mapped reads
  arrange(desc(X3)) %>% 
  head(n = 1) %>% 
  pull(X1)

minor_tmp <- summary$Subtype[2] 
minor_ref <- df %>% 
  filter(Subtype == minor_tmp) %>% 
  # Choose the reference with most mapped reads
  arrange(desc(X3)) %>% 
  head(n = 1) %>% 
  pull(X1)
  
# Write files
write_lines(major_ref, file = paste0(sampleName, ".major_ref.txt"))
write_lines(minor_ref, file = paste0(sampleName, ".minor_ref.txt"))

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")