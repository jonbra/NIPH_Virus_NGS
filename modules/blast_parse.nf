process BLAST_PARSE {

    container 'jonbra/tidyverse_seqinr:1.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/5_blast/", mode:'copy', pattern:'*.{tsv,txt,fa}'

    input:
    tuple val(sampleName), path(scaffolds)
    path references

    output:
    tuple val(sampleName), path('*.{txt,fa}')     , optional:true, emit: subtypes
    tuple val(sampleName), path('*.tsv')     , optional:true, emit: genotypes

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(seqinr)

    scaffolds <- read_tsv("/home/jonr/Prosjekter/learning_nextflow/HCV_results/5_blast/Mix_blast.out",
             col_names = FALSE) %>% 
      # Add a column for the genotype
      separate(X2, into = c("genotype", NA), remove = FALSE) %>% 
      # Count the number of each subject
      add_count(X2) 

    # Which genotypes are present?
    genotypes <- scaffolds %>% 
        distinct(genotype)

    write_tsv(genotypes, paste0("${sampleName}", "_genotypes.tsv"))

    # What is the most common subtype per genotype?
    subtypes <- scaffolds %>% 
        group_by(genotype) %>% 
        slice_max(n, n = 1) %>% 
        distinct(X2)

    # Read the reference fasta file
    fasta <- read.fasta(file = "${references}")

    # Write out the name of each subtype
    for (i in 1:nrow(subtypes)) {
        write_tsv(subtypes[i, 1], file = paste0("${sampleName}_", subtypes\$X2[i], ".txt"), col_names = FALSE)
        write.fasta(sequences = fasta[[subtypes\$X2[i]]], names = subtypes\$X2[i], file.out = paste0(subtypes\$X2[i], ".fa"))
    } 
    """
}
