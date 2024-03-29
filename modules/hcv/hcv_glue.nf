process HCV_GLUE {

    errorStrategy 'ignore'

    publishDir "${params.outdir}/10_glue", mode:'copy', pattern: '*.{json}'

    input:
    tuple val(sampleName), path (bam)

    output:
    path "*.json", optional: true, emit: GLUE_json

    script:
    """
     # Copy bam file to a new filename so they are not present in work directory as links.
    # This is for mounting to the docker image later
    cp ${bam} glue_${bam}
    
    # Pull the latest image
    docker pull cvrbioinformatics/gluetools-mysql:latest
    docker pull cvrbioinformatics/gluetools:latest

    # Start the gluetools-mysql containter
    docker start gluetools-mysql 
    docker exec gluetools-mysql installGlueProject.sh ncbi_hcv_glue

    docker run --rm \
       --name gluetools \
        -v \$PWD:/opt/bams \
        -w /opt/bams \
        --link gluetools-mysql \
        cvrbioinformatics/gluetools:latest gluetools.sh \
         -p cmd-result-format:json \
        -EC \
        -i project hcv module phdrReportingController invoke-function reportBam glue_${bam} 15.0 > ${bam}.json
    """
}      