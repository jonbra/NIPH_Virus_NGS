process SUMMARIZE_MAPPING_HCV {
    
    container 'jonbra/tidyverse_seqinr:2.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/8_summaries/", mode:'copy', pattern:'*.{csv}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'

    input:
    path 'stats/'
    path 'depth/'
    path 'blast/'
    path 'json/'

    output:
    path '*csv', emit: mapping_summary
    path '*{log,sh}'

    script:
    """
    summarize_mapping_hcv.R

    cp .command.log summarize_mapping_command.log
    cp .command.sh summarize_mapping_command.sh
    """

}