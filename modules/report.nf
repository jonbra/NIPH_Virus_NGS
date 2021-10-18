process REPORT {
    tag "$sampleName"
    errorStrategy { task.exitStatus == 143 ? 'retry' : 'ignore' }

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    val(source)

    publishDir "${params.outdir}/1_fastq/log", mode: 'link', pattern:'*.{log,sh}'

    script:
    """
    echo ${sampleName} > report.txt
    echo ${params.outdir} >> report.txt
    echo ${task.cpus} >> report.txt
    echo ${read1} >> report.txt
    echo ${read2} >> report.txt
    echo ${source} >> report.txt
    """
}
