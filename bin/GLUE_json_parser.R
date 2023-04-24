#!/usr/bin/env Rscript

library(tidyverse)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: GLUE_json_parser.R <summary.csv>", call. = FALSE)
}


# Til Kamillas pipeline
# Usage: GLUE_json_parser.R summary.csv

json_files <- list.files(path = "GLUE-rapport_json/",
                         pattern = "json$",
                         full.names = TRUE)

# Create final data file
df_final <- tibble(
  "Sample" = character(), 
  "glecaprevir" = character(), 
  "glecaprevir_mut" = character(), 
  "grazoprevir" = character(), 
  "grazoprevir_mut" = character(), 
  "paritaprevir" = character(), 
  "paritaprevir_mut" = character(), 
  "voxilaprevir" = character(), 
  "voxilaprevir_mut" = character(), 
  "NS34A" = character(),
  "daclatasvir" = character(), 
  "daclatasvir_mut" = character(), 
  "elbasvir" = character(), 
  "elbasvir_mut" = character(), 
  "ledipasvir" = character(), 
  "ledipasvir_mut" = character(), 
  "ombitasvir" = character(), 
  "ombitasvir_mut" = character(), 
  "pibrentasvir" = character(),
  "pibrentasvir_mut" = character(),
  "velpatasvir" = character(), 
  "velpatasvir_mut" = character(), 
  "NS5A" = character(),
  "dasabuvir" = character(),
  "dasabuvir_mut" = character(),
  "sofosbuvir" = character(),
  "sofosbuvir_mut" = character(),
  "NS5B" = character(),
  "HCV project version" = character(),
  "GLUE engine version" = character(),
  "PHE drug resistance extension version"  = character()
)

# List of genes and drugs
genes_drugs <- list("NS3"  = c("glecaprevir", "grazoprevir", "paritaprevir", "voxilaprevir"), 
                    "NS5A" = c("daclatasvir", "elbasvir", "ledipasvir", "ombitasvir", "pibrentasvir", "velpatasvir"), 
                    "NS5B" = c("dasabuvir", "sofosbuvir"))

for (x in 1:length(json_files)) {
  json <- read_json(json_files[x])
  
  # Sample name
  #sample <- str_split_1(basename(json_files[x]), "\\.")[1]
  # Kamillas script
  sample <- unlist(strsplit(basename(json_files[x]), "\\."))[1]
  
  # Få tak i subtype
  subtype <- json[["phdrReport"]][["samReferenceResult"]][["genotypingResult"]][["subtypeCladeCategoryResult"]][["shortRenderedName"]]
  
  # Versjoner:
  projectVersion <- json[["phdrReport"]][["projectVersion"]]
  extensionVersion <- json[["phdrReport"]][["extensionVersion"]]
  engineVersion <- json[["phdrReport"]][["engineVersion"]]
  
  # One row per sample
  # Create a temporary dataframe to populate
  df <- as.data.frame(matrix(nrow = 1, ncol = 31))
  colnames(df) <- c("Sample", 
                    "glecaprevir", 
                    "glecaprevir_mut", 
                    "grazoprevir", 
                    "grazoprevir_mut", 
                    "paritaprevir", 
                    "paritaprevir_mut", 
                    "voxilaprevir", 
                    "voxilaprevir_mut", 
                    "NS34A",
                    "daclatasvir", 
                    "daclatasvir_mut", 
                    "elbasvir", 
                    "elbasvir_mut", 
                    "ledipasvir", 
                    "ledipasvir_mut", 
                    "ombitasvir", 
                    "ombitasvir_mut", 
                    "pibrentasvir",
                    "pibrentasvir_mut",
                    "velpatasvir", 
                    "velpatasvir_mut", 
                    "NS5A",
                    "dasabuvir",
                    "dasabuvir_mut",
                    "sofosbuvir",
                    "sofosbuvir_mut",
                    "NS5B",
                    "HCV project version",
                    "GLUE engine version",
                    "PHE drug resistance extension version")
  
  df$Sample[1] <- sample
  df$`HCV project version` <- projectVersion
  df$`GLUE engine version` <- engineVersion
  df$`PHE drug resistance extension version` <- extensionVersion 

  # Dette er underlisten for Drug Scores. Lengden av denne angir hvor mange drugs som er funnet.
  # Under drugScores så er det en ny liste for hver drug category
  if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]]) > 0) {
    for (i in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]])) {
      # Så er det en ny liste innenfor hver drug score igjen med drug for hver kategori. Denne heter drugAssessemnts
      for (k in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]])) {
        for (l in 1:length(genes_drugs)){
          if (str_detect(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["category"]], names(genes_drugs)[l])) {
            for (m in 1:length(genes_drugs[[l]])) {
              if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == genes_drugs[[l]][m] & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) {
                df[[genes_drugs[[l]][m]]] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
                
                # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i den
                mut <- vector(mode = "character") # Create empty vector to hold mutations
                if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
                  for (n in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]])) {
                    mut <- c(mut, json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[n]][["displayStructure"]])
                  }
                } 
                if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
                  for (n in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]])) {
                    mut <- c(mut, json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[n]][["displayStructure"]])
                  }
                } 
                if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
                  for (n in 1:length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]])) {
                    mut <- c(mut, json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[n]][["displayStructure"]])
                  }
                }
                mut <- paste(mut, collapse = ";")
                df[[paste0(genes_drugs[[l]][m], "_mut")]] <- mut
              }
            }
          }
        }
      }
    }
  }
  
  # Then join mutations per drug category
  tmp <- as_tibble(df)
  
  tmp <- tmp %>% 
    unite("NS34A", c(glecaprevir_mut, grazoprevir_mut, paritaprevir_mut, voxilaprevir_mut), sep = ";", na.rm = TRUE) %>% 
    unite("NS5A", c(daclatasvir_mut, elbasvir_mut, ledipasvir_mut, ombitasvir_mut, pibrentasvir_mut, velpatasvir_mut), sep = ";", na.rm = TRUE) %>% 
    unite("NS5B", c(dasabuvir_mut, sofosbuvir_mut), sep = ";", na.rm = TRUE) 
  
  # Gjøre om innholdet i cellene til en vector, deretter fjerne dupliater i vektoren
  try(df$NS34A <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS34A, ";")), ","))), "\\+"))), collapse = ";"))
  #try(df$NS34A <- unique(unlist(strsplit(tmp %>% pull(NS34A), ";"))))
  try(df$NS5A <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS5A, ";")), ","))), "\\+"))), collapse = ";"))
  #try(df$NS5A <- unique(unlist(strsplit(tmp %>% pull(NS5A), ";"))))
  try(df$NS5B <- paste(unique(unlist(strsplit(gsub(" ", "", unlist(strsplit(unlist(strsplit(tmp$NS5B, ";")), ","))), "\\+"))), collapse = ";"))
  #try(df$NS5B <- unique(unlist(strsplit(tmp %>% pull(NS5B), ";"))))
  
  df <- as_tibble(df)
  
  # Merge with final data structure
  df_final <- bind_rows(df_final, df)
}

# Kamillas script: Må binde df_final til summary-fila igjen
summary <- read_tsv(args[1]) %>% # summary <- read_tsv("/home/jonr/Prosjekter/learning_nextflow/json-test/Run837_HCV_summaries_v7sort.tsv")
  rename("Sample" = "Parameters:")

# join the data
summary_final <- left_join(summary, df_final)

write_tsv(summary_final, file = "summary_with_glue.tsv")

