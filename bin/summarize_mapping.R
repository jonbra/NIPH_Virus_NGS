#!/usr/bin/env Rscript

library(tidyverse)

# Number of mapped reads --------------------------------------------------
path_1 <- "stats/"
path_2 <- "depth/"

df <- list.files(path = path_1, pattern = "stats", full.names = TRUE) %>% 
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
write_csv(df, file = "Genotype_mapping_summary_wide.csv")


# Coverage ----------------------------------------------------------------

# Thoughts: Calculate the number of sites with a coverage at least 5

# List files
files <- list.files(path = path_2, pattern = "gz", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(files), ncol = 3))
colnames(tmp_df) <- c("sampleName", "reference", "cov_breadth")

for (i in 1:length(files)) {
  
  # Get sample name
  tmp_df$sampleName[i] <- str_split(basename(files[i]), "\\.")[[1]][1]
  
  # Get reference name
  tmp_df$reference[i] <- str_remove(str_split(basename(files[i]), "\\.")[[1]][3], "_ref")
  
  # Read the depth per position
  cov <- read_tsv(files[i], col_names = FALSE) 
  
  # Reference length
  ref_length <- nrow(cov)
  
  # Nr. of positions with coverage >= 5
  pos <- nrow(
    cov %>% 
      filter(X3 > 4)
  )
  
  # Coverage breadth
  breadth <- round(pos / ref_length * 100, digits = 2)
  tmp_df$cov_breadth[i] <- breadth
  
}

# Create column for genotype
tmp_df <- as_tibble(tmp_df)
tmp_df <- tmp_df %>% 
  separate(reference, into = c("genotype", NA), sep = "_", remove = FALSE) %>% 
  select(sampleName, genotype, reference, cov_breadth)

# Add the number of mapped reads
# NB: Temporary cleaning of sample names
tmp <- df %>% 
  mutate(sampleName = basename(sampleName)) %>% 
  pivot_longer(!sampleName, names_to = "reference", values_to = "count") %>% filter(!is.na(count))
tmp_df <- left_join(tmp_df, tmp)

# Write file
write_csv(tmp_df, file = "Genotype_mapping_summary_long.csv")

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")

