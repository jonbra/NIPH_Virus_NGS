#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
  stop("Usage: blast_parse.R <sampleName> <blast_out> <scaffolds> <references> <agens>", call.=FALSE)
}

path_1 <- "stats/"
path_2 <- "depth/"
path_3 <- "blast/"
agens  <- args[1]
  
  # Number of mapped reads --------------------------------------------------
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

# List files
files <- list.files(path = path_2, pattern = "gz", full.names = TRUE)

# Empty df
tmp_df <- as.data.frame(matrix(nrow = length(files), ncol = 3))
colnames(tmp_df) <- c("sampleName", "reference", "cov_breadth_min_5")

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
  tmp_df$cov_breadth_min_5[i] <- breadth
  
}

# Create column for genotype
tmp_df <- as_tibble(tmp_df)
tmp_df <- tmp_df %>% 
  separate(reference, into = c("genotype", NA), sep = "_", remove = FALSE) %>% 
  select(sampleName, genotype, reference, cov_breadth_min_5)

# Add the number of mapped reads
# NB: Temporary cleaning of sample names
tmp <- df %>% 
  mutate(sampleName = basename(sampleName)) %>% 
  pivot_longer(!sampleName, names_to = "reference", values_to = "mapped_reads") %>% filter(!is.na(mapped_reads))
tmp_df <- left_join(tmp_df, tmp)

# Create heatmap

if (agens == "ROV") {
  tmp_df %>% 
    # Clean up segment names
    mutate("segment" = case_when(
      genotype == "ASeg1VP1" ~ "Seg1_VP1",
      genotype == "ASeg10NSP4" ~ "Seg10_NSP4",
      genotype == "ASeg2VP2" ~ "Seg2_VP2",
      genotype == "ASeg3VP3" ~ "Seg3_VP3",
      genotype == "ASeg4VP4" ~ "Seg4_VP4_P",
      genotype == "ASeg5NSP1" ~ "Seg5_NSP1",
      genotype == "ASeg6VP6" ~ "Seg6_VP6",
      genotype == "ASeg8NSP2" ~ "Seg8_NSP2",
      genotype == "ASeg9VP7" ~ "Seg9_VP7_G",
      genotype == "BatNSP3" ~ "Seg7_NSP3_Bat",
      genotype == "ASeg11NSP5" ~ "Seg11_NSP5",
      genotype == "ASeg7NSP3" ~ "Seg7_NSP3",
      genotype == "BatVP2" ~ "Seg2_VP2_Bat",
      genotype == "DChickSeg5" ~ "Seg5_NSP1_D"
    )) %>% 
    # Create a factor for ordering
    mutate(segment_f = factor(
      x = segment,
      levels = c("Seg1_VP1", "Seg2_VP2", "Seg2_VP2_Bat", "Seg3_VP3", "Seg4_VP4_P", "Seg5_NSP1", "Seg5_NSP1_D", "Seg6_VP6", "Seg7_NSP3", "Seg7_NSP3_Bat", "Seg8_NSP2", "Seg9_VP7_G", "Seg10_NSP4", "Seg11_NSP5")
    )) %>% 
  ggplot() +
    aes(x = segment_f, y = sampleName) + 
    geom_tile(aes(fill = cov_breadth_min_5)) +
    # Label with coverage and mapped reads
    geom_text(aes(label = paste0(cov_breadth_min_5, "\n", mapped_reads)), color = "white") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
  
  # Write plot
  ggsave("ROV_heatmap_cov_breadth.png", device = "png")
}


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

# Add scaffold lengths to tmp_df
tmp_df <- left_join(tmp_df, df)

# Write file
write_csv(tmp_df, file = "Genotype_mapping_summary_long.csv")

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")

