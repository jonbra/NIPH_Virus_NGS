#!/usr/bin/env Rscript

library(tidyverse)

# Read a list of input files into a list
df <- list.files(path = "plots/", pattern = "gz", full.names = TRUE) %>%
  # Keep the file names as the names of the list elements
  set_names() %>% 
  map(read_tsv, col_names = FALSE) %>% 
  # Reduce the list to a single dataframe. Keep the filenames (list element names) in column 1
  # The column name will be "sampleName"
  bind_rows(.id = "sampleName") %>% 
  # Clean up sampleName
  mutate(sampleName = str_remove(sampleName, "plots//")) %>% 
  # Split on the first "."
  separate(sampleName, into = c("sampleName"), sep = "\\.") %>% 
  # Create a new column that joints the sampleName and reference name
  unite("Plot_name", c("sampleName", "X1"), sep = ".", remove = FALSE) %>% 
  # Rename columns
  rename("Genotype" = X1,
         "Position" = X2,
         "Coverage" = X3)


plots <- df %>% 
  group_by(Plot_name) %>% 
  group_map(
    ~ ggplot(.) + 
      aes(x = Position, y = Coverage) + 
      geom_line() + 
      geom_hline(yintercept = 10, color = "darkgreen", linetype = "dotted") +
      annotate("text", x=900, y=-150, label="Coverage cutoff = 10") +
      ggtitle(.y[[1]]) # .y contains the grouping variable. I.e. the sampleName in this case
    )

# Save the plots
for(i in 1:length(plots)){
  ggsave(plot = plots[[i]], 
         file = paste0((df %>% distinct(Plot_name) %>% pull(Plot_name))[i], # This pulls the ith element of the sampleNames from df
                       i, 
                       ".png"), 
         device = "png", 
         dpi = 300)
}

# Write out sessionInfo() to track versions
session <- capture.output(sessionInfo())
write_lines(session, file = "R_versions.txt")