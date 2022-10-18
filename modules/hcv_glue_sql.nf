process HCV_GLUE_SQL {

    publishDir "${params.outdir}/7_glue", mode:'copy', pattern: '*.{html}'

    input:
    tuple val(sampleName), path(scaffolds)

    output:
    path "*.html"

    script:
    """
    # Pull the latest images
    docker pull cvrbioinformatics/gluetools-mysql:latest
    docker pull cvrbioinformatics/gluetools:latest

    docker exec gluetools-mysql installGlueProject.sh ncbi_hcv_glue
#    docker run --rm \
#       --name gluetools \
#        -v \$PWD/${scaffolds}:/opt/input/scaffolds.fa \
#        -w /opt/input \
#        --link gluetools-mysql \
#        cvrbioinformatics/gluetools:latest gluetools.sh \
#        --console-option log-level:FINEST \
#        --inline-cmd project hcv module phdrReportingController invoke-function reportFastaAsHtml /opt/input/scaffolds.fa ${sampleName}.html

    docker run --rm \
       --name gluetools \
        -v /home/jonr/NODE_1.fa:/opt/input/scaffolds.fa \
        -v \$PWD:/output \
        -w /opt/input \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
        --console-option \
        log-level:FINEST \
        --inline-cmd project hcv module phdrReportingController invoke-function reportFastaAsHtml /opt/input/scaffolds.fa /output/${sampleName}.html
    """
}