process KRAKEN2 {

    container 'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0'

    publishDir "${params.outdir}/3_kraken2/", mode:'copy', pattern:'*.{txt,yml}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'

    label 'large'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path kraken_all, stageAs: 'db'

    output:
    tuple val(sampleName), path("${sampleName}.kraken2_all.report.txt"), emit: report
    path "*.{log,sh}"

    script:
    """
    kraken2 \\
        --db "db" \\
        --threads $task.cpus \\
        --report ${sampleName}.kraken2_all.report.txt \\
        --paired \\
        $read1 $read2

    cp .command.log ${sampleName}.kraken2_all.log
    cp .command.sh ${sampleName}.kraken2_all.sh
    """
}
