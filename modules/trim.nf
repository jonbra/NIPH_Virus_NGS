process TRIM {
  
    conda 'bioconda::cutadapt=4.2'
    container 'quay.io/biocontainers/cutadapt:3.7--py37h8902056_1'

    publishDir "${params.outdir}/2_trimmed", mode:'copy', pattern:'*.{fastq}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*process_command_cutadapt.log'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.sh'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.yml'

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    tuple val(sampleName), path ("${sampleName}_trimmed_R1.fastq"), path ("${sampleName}_trimmed_R2.fastq"), emit: TRIM_out
    tuple val(sampleName), path ("${sampleName}.cutadapt.log")                                             , emit: log
    path "${sampleName}.process_command_cutadapt.log"
    path "*.{sh,yml}"

    script:
    """
    cutadapt -o ${sampleName}_trimmed_R1.fastq -p ${sampleName}_trimmed_R2.fastq -q 30,30 -m 50 ${read1} ${read2} > ${sampleName}.cutadapt.log

    cp .command.log ${sampleName}.process_command_cutadapt.log
    cp .command.sh ${sampleName}.process_command_cutadapt.sh

    cat <<-END_VERSIONS > cutadapt_versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
