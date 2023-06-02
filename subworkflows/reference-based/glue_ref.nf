process GLUE_REF_BASED {

    publishDir "${params.outdir}/ref-based", mode:'copy', pattern: '*.{html,json}'

    input:
    path 'bams/'

    output:
    path "*.html"
    path "*.json", emit: GLUE_json

    script:
    """
    # Copy bam files from bams/ directory so they are not present in work directory as links.
    # This is for mounting to the docker image later
    cp bams/*.bam .
    
    # Pull the latest image
    docker pull cvrbioinformatics/gluetools-mysql:latest
    docker pull cvrbioinformatics/gluetools:latest

    # Start the gluetools-mysql containter
    docker start gluetools-mysql 
    docker exec gluetools-mysql installGlueProject.sh ncbi_hcv_glue

    # Make a for loop over all consensus-sequences

    for bam in \$(ls *.bam)
    do
    docker run --rm \
       --name gluetools \
        -v \$PWD:/opt/bams \
        -w /opt/bams \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
        --console-option log-level:FINEST \
        --inline-cmd project hcv module phdrReportingController invoke-function reportBamAsHtml \${bam} 15.0 \${bam}.html

    docker run --rm \
       --name gluetools \
        -v \$PWD:/opt/bams \
        -w /opt/bams \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
        --console-option cmd-result-format:json -EC \
        --inline-cmd project hcv module phdrReportingController invoke-function reportBam \${bam} 15.0 > \${bam}.json
    done
    """
}      