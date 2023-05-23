process FASTQC {

    conda "bioconda::fastqc=0.11.9"
    container 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'

    publishDir "${params.outdir}/1_fastqc/", mode:'copy', pattern:'*.{zip}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.yml'

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    val(source)

    output:
    tuple val(sampleName), path ("*.zip"), emit: FASTQC_out
    path "*.{log,sh,yml}"

    script:
    """
    fastqc -t $task.cpus ${read1} ${read2}
    
    cp .command.log ${sampleName}.${source}.fastqc.log
    cp .command.sh ${sampleName}.${source}.fastqc.sh

    cat <<-END_VERSIONS > fastqc_versions.yml
        "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}
