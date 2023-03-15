#!/usr/bin/env Rscript

library(tidyverse)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: GLUE_json_parser.R <summary.csv>", call. = FALSE)
}


# Til Kamillas pipeline
# Usage: GLUE_json_parser.R summary.csv

json_files <- list.files(pattern = "json$",
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
  
# Loop through all the json files and extract relevant info
for (x in 1:length(json_files)) {
  json <- read_json(json_files[x])
  
  # Sample name
  #sample <- str_split_1(basename(json_files[x]), "\\.")[1]
  # Kamillas script
  sample <- str_split_1(basename(json_files[x]), "_")[1]
  
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
        # NS34A category
        if (str_detect(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["category"]], "NS3")) {
          
          # glecaprevrir drug & length > 0 hvis det er insufficient coverage
          if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "glecaprevir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) { 
            df$glecaprevir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i den
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$glecaprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$glecaprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$glecaprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
            # grazoprevir drug  
          } else if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "grazoprevir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) {
            df$grazoprevir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i denne kategorien
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[1]][["rasScores_category_I"]]) > 0) {
              df$grazoprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$grazoprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]] 
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$grazoprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
            # paritaprevir drug 
          } else if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "paritaprevir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) {
            df$paritaprevir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i denne kategorien
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$paritaprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$paritaprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$paritaprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
            # voxilaprevir drug 
          } else if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "voxilaprevir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) {
            df$voxilaprevir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i denne kategorien
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$voxilaprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$voxilaprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$voxilaprevir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
          }
          
          # NS5A category
        } else if (str_detect(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["category"]], "NS5A")) {
          
          # daclatasvir drug
          if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "daclatasvir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) { 
            df$daclatasvir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i den
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$daclatasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$daclatasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$daclatasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
            # elbasvir drug  
          } else if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "elbasvir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) {
            df$elbasvir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i denne kategorien
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$elbasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$elbasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$elbasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
            # ledipasvir drug 
          } else if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "ledipasvir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) {
            df$ledipasvir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i denne kategorien
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$ledipasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$ledipasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$ledipasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
            # ombitasvir drug 
          } else if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "ombitasvir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) {
            df$ombitasvir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i denne kategorien
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$ombitasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$ombitasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$ombitasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
            # pibrentasvir drug 
          } else if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "pibrentasvir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) {
            df$pibrentasvir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i denne kategorien
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$pibrentasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$pibrentasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$pibrentasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
            # velpatasvir drug 
          } else if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "velpatasvir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) {
            df$velpatasvir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i denne kategorien
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$velpatasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$velpatasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$velpatasvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
          } 
          
          # NS5B category
        } else if (str_detect(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["category"]], "NS5B")) {
          
          # dasabuvir drug
          if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "dasabuvir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) { 
            df$dasabuvir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i den
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$dasabuvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$dasabuvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$dasabuvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
            }
            
            # sofosbuvir drug  
          } else if (json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drug"]][["id"]] == "sofosbuvir" & length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]) > 0) {
            df$sofosbuvir[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["drugScoreDisplayShort"]]
            
            # De tre kategoriene er lister. Hvis lengden er > 0 betyr det at det er en mutasjon i denne kategorien
            if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]]) > 0) {
              df$sofosbuvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_I"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]]) > 0) {
              df$sofosbuvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_II"]][[1]][["displayStructure"]]
            } else if (length(json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]]) > 0) {
              df$sofosbuvir_mut[1] <- json[["phdrReport"]][["samReferenceResult"]][["drugScores"]][[i]][["drugAssessments"]][[k]][["rasScores_category_III"]][[1]][["displayStructure"]]
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
  df$NS34A <- unique(str_split_1(tmp %>% pull(NS34A), ";"))
  df$NS5A <- unique(str_split_1(tmp %>% pull(NS5A), ";"))
  df$NS5B <- unique(str_split_1(tmp %>% pull(NS5B), ";"))
  
  df <- as_tibble(df)
  
  # Merge with final data structure
  df_final <- bind_rows(df_final, df)
  
}

# Kamillas script: Må binde df_final til summary-fila igjen
sumamry <- read_tsv(args[1]) %>% # summary <- read_tsv("/home/jonr/Prosjekter/learning_nextflow/json-test/Run837_HCV_summaries_v7sort.tsv")
  rename("Sample" = "Parameters:")

# join the data
summary_final <- left_join(summary, df_final)

write_tsv(summary_final, file = "summary_with_glue.tsv")

