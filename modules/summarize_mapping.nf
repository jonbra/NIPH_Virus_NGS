process SUMMARIZE_MAPPING {
    
    container 'jonbra/tidyverse_seqinr:1.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/8_summaries/", mode:'copy', pattern:'*.{csv}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.txt'

    input:
    path 'stats/'

    output:
    path '*csv', emit: mapping_summary
    path 'R_versions.txt'
    path '*{log,sh}'

    script:
    """
    summarize_mapping.R

    cp .command.log summarize_mapping_command.log
    cp .command.sh summarize_mapping_command.sh
    """

}