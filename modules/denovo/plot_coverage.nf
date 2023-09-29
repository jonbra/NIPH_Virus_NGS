process PLOT_COVERAGE {
    
    container 'jonbra/tidyverse_seqinr:2.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/7_coverage/", mode:'copy', pattern:'*.{png}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'

    input:
    path 'plots/'

    output:
    path '*png', emit: coverage_plots
    path '*{log,sh}'
    path "bam_coverage_versions.yml", emit: versions

    script:
    """
    bam_coverage.R

    cp .command.log plot_coverage_command.log
    cp .command.sh plot_coverage_command.sh

    cat <<-END_VERSIONS > bam_coverage_versions.yml
    "${task.process}":
        R-session: \$(cat R_versions.txt)
    END_VERSIONS
    """

}