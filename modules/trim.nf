process TRIM {
  
    container 'quay.io/biocontainers/cutadapt:3.7--py37h8902056_1'

    publishDir "${params.outdir}/2_trimmed", mode:'copy', pattern:'*.{fastq}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.yml'

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    tuple val(sampleName), path ("${sampleName}_trimmed_R1.fastq"), path ("${sampleName}_trimmed_R2.fastq"), emit: TRIM_out
    path "*.{log,sh,yml}"

    script:
    """
    cutadapt -o ${sampleName}_trimmed_R1.fastq -p ${sampleName}_trimmed_R2.fastq -q 30,30 -m 50 ${read1} ${read2}

    cp .command.log ${sampleName}.cutadapt.log
    cp .command.sh ${sampleName}.cutadapt.sh

    cat <<-END_VERSIONS > cutadapt_versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
