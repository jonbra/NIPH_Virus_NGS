process BLAST_PARSE {

    container 'jonbra/tidyverse_seqinr:1.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/5_blast/", mode:'copy', pattern:'*.{tsv,txt,fa}'

    input:
    //tuple val(sampleName), path(blast_out)
    path 'blast/*'
    path references
    //tuple val(sampleName), path(scaffolds_fasta)
    path 'scaffolds/*'

    output:
    tuple val(sampleName), path('*.{txt}'), emit: subtypes
    tuple val(sampleName), path('*.tsv')  , emit: genotypes
    path '*ref.fa'                        , emit: subtypes_references
    path '*scaffolds.fa'                  , emit: scaffolds_fasta

    script:
    """
    NB! Det blir ulike sample names in!!
    #!/usr/bin/env Rscript

    library(tidyverse)
    library(seqinr)

    scaffolds <- read_tsv("${blast_out}",
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
        write.fasta(sequences = fasta[[subtypes\$X2[i]]], names = subtypes\$X2[i], file.out = paste0("${sampleName}_", subtypes\$X2[i], "_ref.fa"))
    } 

    # NB! Dette er spades scaffolds som blir lest?? Blir ikke riktig å splitte slik?
    # Split scaffolds per subtype
    # Read scaffolds fasta
    scaffolds_fa <- read.fasta(file = "${scaffolds_fasta}")

    # Split scaffolds per genotype
    split <- scaffolds %>% 
        group_split(genotype)

    # Write one fasta per genotype
    for (i in 1:length(split)){
        # Store the scaffold names and corresponding genotype
        tmp <- split[[i]] %>% select(X1, genotype)
        # Subset the scaffolds
        geno_fa <- scaffolds_fa[tmp\$X1]
        # Write to file
        write.fasta(sequences = geno_fa, names = names(geno_fa), file.out = paste0("${sampleName}_", tmp\$genotype[1], "_scaffolds.fa"))
    }

    # Write out sessionInfo() to track versions
    session <- capture.output(sessionInfo())
    write_lines(session, file = "R_versions.txt")
    """
}
