process MULTIQC {

    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'

    label 'small'

    input:
    path 'data/*'

    output:
    path '0_multiqc_report.html'
    path '0_multiqc_data'
    path 'all.multiqc.{log,sh}'

    publishDir "${params.outdir}/6_multiqc/", mode: 'copy', pattern:'*multiqc_*'
    publishDir "${params.outdir}/6_multiqc/", mode: 'copy', pattern:'*.yml'

    script:
    """
    multiqc data
    mv multiqc_data 0_multiqc_data
    mv multiqc_report.html 0_multiqc_report.html

    cp .command.sh all.multiqc.sh
    cp .command.log all.multiqc.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
