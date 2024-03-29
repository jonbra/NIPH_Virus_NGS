process BLAST_PARSE {

    // Create a custom conda environment with the necessary R packages?
    container 'jonbra/tidyverse_seqinr:2.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/5_blast/", mode:'copy', pattern:'*.{csv,tsv,fa}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'

    input:
    tuple val(sampleName), path(blast_out), path(scaffolds), path(read1), path(read2), path(references)
    //val references

    output:
    tuple val(sampleName), path("${sampleName}*ref.fa"), path(read1), path(read2)          , emit: FOR_MAPPING
    tuple val(sampleName), path("${sampleName}*ref.fa"), path("${sampleName}*scaffolds.fa"), emit: FOR_ABACAS
    tuple val(sampleName), path('*.csv')                                                   , emit: blast_res
    path '*scaffolds.fa'                                                                   , emit: RESISTANCE_BLAST
    path "*.{log,sh}"
    path "blast_parse_versions.yml", emit: versions

    script:
    """
    # params.agens comes from the agent-specific config file
    blast_parse.R "$sampleName" "$blast_out" "$scaffolds" "$references" $params.agens

    cp .command.log ${sampleName}.blast_parse_command.log
    cp .command.sh ${sampleName}.blast_parse_command.sh

    cat <<-END_VERSIONS > blast_parse_versions.yml
    "${task.process}":
        R-session: \$(cat R_versions.txt)
    END_VERSIONS
    """
}
