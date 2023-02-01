process HCV_GLUE_SQL {

    container 'docker:latest'

    publishDir "${params.outdir}/7_glue", mode:'copy', pattern: '*.{html}'

    input:
    tuple val(sampleName), path("${sampleName}*.bam"), path("${sampleName}*.bai")

    output:
    path "*.html"

    script:
    """
    # Pull the latest images
    docker pull cvrbioinformatics/gluetools-mysql:latest
    docker pull cvrbioinformatics/gluetools:latest

    # Start the gluetools-mysql containter
    docker exec gluetools-mysql installGlueProject.sh ncbi_hcv_glue

    # Make a for loop over all consensus-sequences

    for bam in \$(ls *.bam)
    do
    docker run --rm \
       --name gluetools \
        -v \$PWD/:/opt/input/\${bam} \
        -v \$PWD:/output \
        -w /opt/input \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
        --console-option \
        log-level:FINEST \
        --inline-cmd project hcv module phdrReportingController invoke-function reportBamAsHtml /opt/input/\${bam} /output/\${bam}.html
    done
    """
}      