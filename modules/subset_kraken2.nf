process SUBSET_KRAKEN2 {

    container 'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0'

    publishDir "${params.outdir}/3_kraken2/", mode:'copy', pattern:'*.{fastq}'

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path report
    path kraken_out

    output:
    tuple val(sampleName), path('*taxid*R1.fastq*'), path('*taxid*R2.fastq*'), optional:true, emit: subset_reads_fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    extract_kraken_reads.py \\
        -k $kraken_out \\
        -r $report \\
        -t ${params.taxid} \\
        --include-children \\
        -s1 $read1 \\
        -s2 $read2 \\
        -o ${sampleName}_taxid_${params.taxid}_R1.fastq \\
        -o2 ${sampleName}_taxid_${params.taxid}_R2.fastq \\
        --fastq-output
    """

}