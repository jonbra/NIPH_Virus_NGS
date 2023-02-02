process HCV_GLUE_SQL {

    publishDir "${params.outdir}/7_glue", mode:'copy', pattern: '*.html'

    input:
    path 'bams/'

    output:
    path "*.html"

    script:
    """
    # Copy bam files from bams/ directory so they are not present in work directory as links.
    # This is for mounting to the docker image later
    cp bams/*nodups.bam .
    
    # Pull the latest image
    docker pull cvrbioinformatics/gluetools-mysql:latest
    docker pull cvrbioinformatics/gluetools:latest

    # Start the gluetools-mysql containter
    docker start gluetools-mysql 
    docker exec gluetools-mysql installGlueProject.sh ncbi_hcv_glue

    # Make a for loop over all consensus-sequences

    for bam in \$(ls *nodups.bam)
    do
    docker run --rm \
       --name gluetools \
        -v \$PWD:/opt/bams \
        -w /opt/bams \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
        --console-option log-level:FINEST \
        --inline-cmd project hcv module phdrReportingController invoke-function reportBamAsHtml \${bam} 15.0 \${bam}.html
    done
    """
}      