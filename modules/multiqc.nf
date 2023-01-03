process MULTIQC {

    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'

    label 'small'

    input:
    path 'data/*'

    output:
    path '0_multiqc_report.html'
    path '0_multiqc_data'
    path 'multiqc_command.{log,sh}'
    path 'multiqc_versions.yml'

    publishDir "${params.outdir}/6_multiqc/", mode: 'copy', pattern:'*multiqc_*'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.yml'

    script:
    """
    multiqc data
    mv multiqc_data 0_multiqc_data
    mv multiqc_report.html 0_multiqc_report.html

    cp .command.sh multiqc_command.sh
    cp .command.log multiqc_command.log

    cat <<-END_VERSIONS > multiqc_versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
