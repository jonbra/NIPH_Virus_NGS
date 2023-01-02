library(tidyverse)

# Read Opportuni-C samples
resultater <- readxl::read_xlsx("/mnt/N/Virologi/NGS/1-NGS-Analyser/3-Prosjekt/HCV-prosjekter/OpportuniC/2_Resultater_HCV_OpportuniC.xlsx",
                                sheet = "Resultater",
                                skip = 1)

resultater <- readxl::read_xlsx("/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/HCV/2022/Resultater_HCV_oversikt_2022.xlsx",
                                sheet = "Resultater",
                                skip = 1)
resultater <- resultater %>% 
  rename("Sample_ID" = `Sample_ID\r\n(markert gult de som ikke har full resistensrapport)`)

# List all fastq files
fastq_files_Run613 <- list.files("/mnt/N/NGS/3-Sekvenseringsbiblioteker/2021/Illumina_RunXXX/Run613_Corona_HCV/",
                               recursive = TRUE,
                               full.names = TRUE,
                               pattern = "gz$")
fastq_files_Run622 <- list.files("/mnt/N/NGS/3-Sekvenseringsbiblioteker/2021/Illumina_RunXXX/Run622/",
                                 recursive = TRUE,
                                 full.names = TRUE,
                                 pattern = "gz$")
fastq_files_Run630 <- list.files("/mnt/N/NGS/3-Sekvenseringsbiblioteker/2021/Illumina_RunXXX/Run630/",
                                 recursive = TRUE,
                                 full.names = TRUE,
                                 pattern = "gz$")
fastq_files_Run652 <- list.files("/mnt/N/NGS/3-Sekvenseringsbiblioteker/2021/Illumina_RunXXX/Run652/",
                                 recursive = TRUE,
                                 full.names = TRUE,
                                 pattern = "gz$")
fastq_files_Run712 <- list.files("/mnt/N/NGS/3-Sekvenseringsbiblioteker/2022/Illumina_RunXXX/Run712_Corona+HCV/",
                                 recursive = TRUE,
                                 full.names = TRUE,
                                 pattern = "gz$")
fastq_files_Run726 <- list.files("/mnt/N/NGS/3-Sekvenseringsbiblioteker/2022/Illumina_RunXXX/Run726_Corona_HCV_test_nye_indexer/",
                                 recursive = TRUE,
                                 full.names = TRUE,
                                 pattern = "gz$")
fastq_files_Run738 <- list.files("/mnt/N/NGS/3-Sekvenseringsbiblioteker/2022/Illumina_RunXXX/Run738_Corona_HCV/",
                                 recursive = TRUE,
                                 full.names = TRUE,
                                 pattern = "gz$")
fastq_files_Run820 <- list.files("/mnt/N/NGS/3-Sekvenseringsbiblioteker/2022/Illumina_RunXXX/Run820_Virus/",
                                 recursive = TRUE,
                                 full.names = TRUE,
                                 pattern = "gz$")
fastq_files_all <- c(fastq_files_Run613, fastq_files_Run622, fastq_files_Run630, fastq_files_Run652, fastq_files_Run712, fastq_files_Run726, fastq_files_Run738)
fastq_files_all <- fastq_files_Run820

R1 <- sort(fastq_files_all[grep("R1", fastq_files_all)])
R2 <- sort(fastq_files_all[grep("R2", fastq_files_all)])

df <- as_tibble(cbind(R1, R2))

df <- df %>%
  mutate(sample_id = gsub("_.*", "", basename(R1))) %>%
  select("sample" = sample_id,
         "fastq_1" = R1,
         "fastq_2" = R2)

# Match with the result sheet 
resultater_final <- left_join(resultater, df, by = c("Sample_ID" = "sample")) %>% 
  drop_na(Sample_ID) %>% 
  select("sample" = Sample_ID, fastq_1, fastq_2)

write_csv(resultater_final, "/home/jonr/Prosjekter/learning_nextflow/2022.11.30-HCV_Run820_samplesheet.csv")
