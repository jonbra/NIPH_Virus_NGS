process HCV_KRAKEN2 {

    // This script classifies the input reads with Kraken2 against a custom HCV database.
    // The database was build using the file HCVgenosubtypes_8.5.19_clean.fa which is the same as used in the Blastn
    publishDir "${params.outdir}/hcv_kraken2/", mode:'copy', pattern:'*.{fastq,txt,yml}'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path  kraken_hcv_db

    output:
    tuple val(sampleName), path('*hcv_classified_1*'), path('*hcv_classified_2*')     , optional:true, emit: classified_reads_fastq
    path "*hcv_classifiedreads*", emit: classified_reads_assignment
    path "*hcv_report.txt", emit: hcv_report
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    kraken2 \\
        --db $kraken_hcv_db \\
        --threads $task.cpus \\
        --report ${sampleName}.kraken2.hcv_report.txt \\
        --paired \\
        --classified-out ${sampleName}.hcv_classified#.fastq \\
        --output ${sampleName}.kraken2.hcv_classifiedreads.txt \\
        $read1 $read2


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}