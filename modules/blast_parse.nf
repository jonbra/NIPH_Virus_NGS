process BLAST_PARSE {

    container 'jonbra/tidyverse_seqinr:1.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/5_blast/", mode:'copy', pattern:'*.{tsv,txt,fa}'

    input:
    tuple val(sampleName), path(blast_out), path(scaffolds), path(read1), path(read2)
    path references

    output:
    tuple val(sampleName), path('*ref.fa'), path(read1), path(read2), emit: FOR_MAPPING
    //tuple val(sampleName), path('*.txt')  , emit: subtypes
    tuple val(sampleName), path('*.tsv')  , emit: genotypes
    path '*ref.fa'                        , emit: subtypes_references
    path '*scaffolds.fa'                  , emit: scaffolds_fasta
    path 'R_versions.txt'

    script:
    """
    blast_parse.R "$sampleName" "$blast_out" "$scaffolds" "$references"
    """
}
