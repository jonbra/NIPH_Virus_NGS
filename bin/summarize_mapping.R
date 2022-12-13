#!/usr/bin/env Rscript

library(tidyverse)

df <- list.files(path = "stats/", pattern = "stats", full.names = TRUE) %>% 
  # Keep the file names as the names of the list elements
  set_names() %>% 
  map(read_tsv, col_names = FALSE) %>% 
  # Reduce the list to a single dataframe. Keep the filenames (list element names) in column 1
  # The column name will be "sampleName"
  bind_rows(.id = "sampleName") %>% 
  # Clean up sampleName
  mutate(sampleName = str_remove(sampleName, "stats//")) %>% 
  separate(sampleName, into = c("sampleName", NA, "reference"), sep = "\\.") %>% 
  mutate(reference = str_remove(reference, "_ref")) %>% 
  # Create a new column that joints the sampleName and reference name
  unite("Sample_ref", c("sampleName", "reference"), sep = ".", remove = FALSE) %>% 
  # Rename columns
  rename("Info" = X1,
         "Number" = X2,
         "Comment" = X3) %>% 
  # Extract relevant info
  filter(Info == "reads mapped:") %>% 
  # Drop the Comment column
  select(-Comment) %>% 
  rename("Reads_mapped" = "Number") %>% 
  select(-Info) %>% 
  # Get each sample on a single line
  pivot_wider(names_from = reference, values_from = Reads_mapped, id_cols = sampleName)

# Write file
write_csv(df, file = "Genotype_mapping_summary.csv")

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")

