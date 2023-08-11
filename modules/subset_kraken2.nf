process SUBSET_KRAKEN2 {

    maxForks 1

    container 'quay.io/biocontainers/krakentools:1.2--pyh5e36f6f_0'

    publishDir "${params.outdir}/3_kraken2/", mode:'copy', pattern:'*.fastq'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path report
    path kraken_out

    output:
    tuple val(sampleName), path('*taxid*R1.fastq*'), path('*taxid*R2.fastq*'), emit: subset_reads_fastq, optional: true
    path "versions.yml"                                                      , emit: versions

    script:
    def VERSION = '1.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
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
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extract_kraken_reads.py: ${VERSION}
    END_VERSIONS
    """

}