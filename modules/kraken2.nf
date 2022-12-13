process KRAKEN2 {

    container 'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0'

    publishDir "${params.outdir}/3_kraken2/", mode:'copy', pattern:'*.{txt,yml}'

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path  kraken_main

    output:
    tuple val(sampleName), path("*report.txt"), emit: report
    path "versions.yml"                       , emit: versions

    script:
    """
    kraken2 \\
        --db $kraken_main \\
        --threads $task.cpus \\
        --report ${sampleName}.kraken2.report.txt \\
        --paired \\
        $read1 $read2


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}