#!/usr/bin/env Rscript

library(tidyverse)

# Number of mapped reads --------------------------------------------------
path_1 <- "stats/"
path_2 <- "depth/"
path_3 <- "blast/"
path_4 <- "json/"

# Reads mapped ------------------------------------------------------------
stats <- list.files(path = path_1, pattern = "\\.stats$", full.names = TRUE) %>% 
  # Keep the file names as the names of the list elements
  set_names() %>% 
  map(read_tsv, col_names = FALSE, comment = "#") %>% 
  # Reduce the list to a single dataframe. Keep the filenames (list element names) in column 1
  # The column name will be "sampleName"
  bind_rows(.id = "sampleName") %>% 
  # Clean up sampleName
  mutate(sampleName = str_remove(sampleName, "stats//")) %>%
  # Extract relevant info
  filter(X2 == "reads mapped:") %>%
  # Create a new column that keeps the sample name, reference for mapping and major/minor
  mutate("Sample_ref" = str_remove(sampleName, "\\.markdup\\.bam\\.stats")) %>% 
  # Keep the reference in a separate column
  separate(sampleName, into = c("sampleName", "reference", "major_minor"), sep = "\\.") %>% 
  # Select relevant columns and rename
  select(sampleName,
         reference,
         major_minor,
         Sample_ref,
         "Reads_mapped" = X3) 
  
# Coverage ----------------------------------------------------------------

# List files
cov_files <- list.files(path = path_2, pattern = "txt\\.gz$", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(cov_files), ncol = 4))
colnames(tmp_df) <- c("sampleName", "reference", "cov_breadth_min_5", "major_minor")

for (i in 1:length(cov_files)) {
  
  # Get sample name
  tmp_df$sampleName[i] <- str_split(basename(cov_files[i]), "\\.")[[1]][1]
  
  # Get reference name
  tmp_df$reference[i] <- str_split(basename(cov_files[i]), "\\.")[[1]][2]
  
  # Get major or minor
  tmp_df$major_minor[i] <- str_split(basename(cov_files[i]), "\\.")[[1]][3]
  
  # Read the depth per position
  cov <- read_tsv(cov_files[i], col_names = FALSE) 
  
  # Reference length
  ref_length <- nrow(cov)
  
  # Nr. of positions with coverage >= 5
  pos <- nrow(
    cov %>% 
      filter(X3 >= 5)
  )
  
  # Coverage breadth
  breadth <- round(pos / ref_length * 100, digits = 2)
  tmp_df$cov_breadth_min_5[i] <- breadth
  
}

# Create column for subtype and Sample_ref
tmp_df <- as_tibble(tmp_df)
tmp_df <- tmp_df %>% 
  separate(reference, into = c("subtype", NA), sep = "_", remove = FALSE) %>%
  unite("Sample_ref", c("sampleName", "reference", "major_minor"), sep = ".", remove = FALSE) %>% 
  select(sampleName, reference, major_minor, Sample_ref, subtype, cov_breadth_min_5)

# Add the number of mapped reads. for the first mapping, the majority mapping (with dups and without dups)
tmp_df <- tmp_df %>%
  # Join if sampleName, reference and major_minor are the same. Get reads mapped without duplicates
  left_join(stats %>% select(-Sample_ref), by = join_by(sampleName, reference, major_minor)) %>% 
  rename("Reads_mapped_no_duplicates" = Reads_mapped)

# Add reads mapped from the first mapping against all references
tmp_df <- stats %>% 
  filter(reference == "first_mapping") %>% 
  select(-major_minor, -Sample_ref, -reference) %>% 
  right_join(tmp_df) %>% 
  rename("Reads_mapped_first_mapping" = Reads_mapped)

# Add reads mapped from second mapping, including duplicates
tmp_df <- stats %>% 
  filter(reference == "majority_mapping") %>% 
  select(-major_minor, -Sample_ref, -reference) %>% 
  right_join(tmp_df) %>% 
  rename("Reads_mapped_second_mapping" = Reads_mapped)

# Length of scaffolds -----------------------------------------------------

# Print the length of the longest scaffold matching the given reference
blast_files <- list.files(path = path_3, pattern = "out$", full.names = TRUE)

# Empty df
df <- tribble(
  ~"scaffold_length", ~"reference", ~"sampleName",
  )

for (i in 1:length(blast_files)) {
  
  # Read the length of scaffolds
  blast_out <- read_tsv(blast_files[i], col_names = FALSE) %>% 
    # Separate the genotype from the subject header
    separate(X2, into = c("genotype", NA), remove = FALSE) %>% 
    # Get scaffold length info
    separate(X1, c(NA, NA, NA, "scaffold_length", NA, NA), sep = "_", remove = FALSE) %>% 
    mutate(scaffold_length = as.numeric(scaffold_length)) %>% 
    # Select the row with the longest scaffold lengths for each genotype/blast hit
    group_by(genotype) %>% 
    slice_max(order_by = scaffold_length, n = 1) %>% 
    ungroup %>% 
    # Sometimes the same scaffolds has two or more hits
    distinct(X1, .keep_all = TRUE) %>% 
    select(scaffold_length,
           "reference" = X2) %>% 
    # Legg til en kolonne med Sample Name
    add_column("sampleName" = str_split(basename(blast_files[i]), "_")[[1]][1])
  
  # Add to df
  df <- bind_rows(df, blast_out)
  
}

# Add scaffold lengths to tmp_df. Both sampleName and reference must match in order to join
tmp_df <- left_join(tmp_df, df)


# GLUE --------------------------------------------------------------------

glue_reports <- list.files(path = path_4, pattern = "GLUE_report.tsv$", full.names = TRUE) %>% 
  # Keep the file names as the names of the list elements
  set_names() %>% 
  map(read_tsv, col_types = cols(GLUE_subtype = col_character())) %>% 
  # Reduce the list to a single dataframe. Keep the filenames (list element names) in column 1
  # The column name will be "sampleName"
  bind_rows(.id = "sampleName") %>% 
  # Clean up sampleName
  mutate(sampleName = str_remove(sampleName, "json//")) %>% # "json//
  # Create a new column that keeps the sample name, reference for mapping and major/minor
  mutate("Sample_ref" = str_remove(sampleName, "_GLUE_report\\.tsv")) %>% 
  # Keep the reference in a separate column
  separate(sampleName, into = c("sampleName", "reference", "major_minor"), sep = "\\.") %>% 
  mutate(major_minor = str_remove(major_minor, "_GLUE_report")) 

# Join on sampleName, reference, major_minor
tmp_df <- left_join(tmp_df, glue_reports)

# Write file
write_csv(tmp_df, file = "Genotype_mapping_summary_long.csv")

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")
