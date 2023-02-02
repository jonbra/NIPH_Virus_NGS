process BLAST_PARSE {

    container 'jonbra/tidyverse_seqinr:1.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/5_blast/", mode:'copy', pattern:'*.{csv,tsv,fa}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.txt'

    input:
    tuple val(sampleName), path(blast_out), path(scaffolds), path(read1), path(read2), path(references)
    //val references

    output:
    tuple val(sampleName), path("${sampleName}*ref.fa"), path(read1), path(read2)          , emit: FOR_MAPPING
    tuple val(sampleName), path("${sampleName}*ref.fa"), path("${sampleName}*scaffolds.fa"), emit: FOR_ABACAS
    tuple val(sampleName), path('*.csv')                                                   , emit: blast_res
    path '*scaffolds.fa'                                                                   , emit: RESISTANCE_BLAST
    path 'R_versions.txt'
    path "*.{log,sh}"

    script:
    """
    # params.agens comes from the agent-specific config file
    blast_parse.R "$sampleName" "$blast_out" "$scaffolds" "$references" $params.agens

    cp .command.log ${sampleName}.blast_parse_command.log
    cp .command.sh ${sampleName}.blast_parse_command.sh
    """
}
