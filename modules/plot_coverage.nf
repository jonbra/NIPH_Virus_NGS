process PLOT_COVERAGE {
    
    container 'jonbra/tidyverse_seqinr:2.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/7_coverage/", mode:'copy', pattern:'*.{png}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.txt'

    input:
    path 'plots/'

    output:
    path '*png', emit: coverage_plots
    path 'R_versions.txt'
    path '*{log,sh}'

    script:
    """
    bam_coverage.R

    cp .command.log plot_coverage_command.log
    cp .command.sh plot_coverage_command.sh
    """

}