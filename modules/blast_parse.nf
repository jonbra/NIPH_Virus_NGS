process BLAST_PARSE {

    container 'rocker/tidyverse:4.2.1'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/5_blast/", mode:'copy', pattern:'*.{txt,tsv}'

    input:
    tuple val(sampleName), path(scaffolds)

    output:
    tuple val(sampleName), path('*.txt')     , optional:true, emit: subtypes
    tuple val(sampleName), path('*.tsv')     , optional:true, emit: genotypes

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    scaffolds <- read_tsv("${scaffolds}",
         col_names = FALSE) %>% 
    # Add a column for the genotype
    separate(X2, into = c("genotype", NA), remove = FALSE) %>% 
    # Count the number of each query
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

    # Write out the name of each subtype
    for (i in 1:nrow(subtypes)) {
    write_tsv(subtypes[i, 1], file = paste0("${sampleName}_", subtypes[i, 1], ".txt"), col_names = FALSE)
    }
 
    """
}
