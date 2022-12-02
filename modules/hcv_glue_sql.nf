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

    # Start the gluetools-mysql containter
    docker exec gluetools-mysql installGlueProject.sh ncbi_hcv_glue

    # Make a for loop over all consensus-sequences

    # Start the genotyping and resistance analysis
    # NB: Not working if multiple scaffolds in the scaffolds file
    docker run --rm \
       --name gluetools \
        -v \$PWD/${scaffolds}:/opt/input/${sampleName}.fa \
        -v \$PWD:/output \
        -w /opt/input \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
        --console-option \
        log-level:FINEST \
        --inline-cmd project hcv module phdrReportingController invoke-function reportFastaAsHtml /opt/input/${sampleName}.fa /output/${sampleName}.html
    """
}      