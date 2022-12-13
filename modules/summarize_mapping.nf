process SUMMARIZE_MAPPING {
    
    container 'jonbra/tidyverse_seqinr:1.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/8_summaries/", mode:'copy', pattern:'*.{csv,txt}'

    input:
    path 'stats/'

    output:
    path '*csv', emit: mapping_summary
    path 'R_versions.txt'

    script:
    """
    summarize_mapping.R
    """

}