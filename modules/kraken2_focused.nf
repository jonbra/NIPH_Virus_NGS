process KRAKEN2_FOCUSED {

    conda "bioconda::kraken2=2.1.2 conda-forge::pigz=2.6"
    container 'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0'

    publishDir "${params.outdir}/3_kraken2/", mode:'copy', pattern:'*.{txt}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.yml'

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path kraken_db

    output:
    tuple val(sampleName), path("${sampleName}.classified_1.fastq"), path("${sampleName}.classified_2.fastq"), emit: classified_reads_fastq
    tuple val(sampleName), path("${sampleName}.kraken2_focused.report.txt")                                  , emit: report
    path "*.{log,sh,yml}"

    script:
    """
    kraken2 \\
        --db $kraken_db \\
        --threads $task.cpus \\
        --report ${sampleName}.kraken2_focused.report.txt \\
        --paired \\
        --classified-out ${sampleName}.classified#.fastq \\
        --output ${sampleName}.kraken2_focused.classifiedreads.txt \\
        $read1 $read2

    cp .command.log ${sampleName}.kraken2_focused.log
    cp .command.sh ${sampleName}.kraken2_focused.sh

    cat <<-END_VERSIONS > kraken2_versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
