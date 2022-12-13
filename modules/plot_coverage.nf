process PLOT_COVERAGE {
    
    container 'jonbra/tidyverse_seqinr:1.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/7_coverage/", mode:'copy', pattern:'*.{png,txt}'

    input:
    path 'plots/'

    output:
    path '*png', emit: coverage_plots
    path 'R_versions.txt'

    script:
    """
    bam_coverage.R
    """

}