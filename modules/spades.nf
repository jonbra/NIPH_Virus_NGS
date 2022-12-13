process SPADES {

    container 'quay.io/biocontainers/spades:3.15.4--h95f258a_0'

    // Sometimes there can be zero reads after SUBSET_KRAKEN2, in which case Spades will crash.
    errorStrategy 'ignore'

    label 'spades'

    publishDir "${params.outdir}/4_spades/", mode:'copy', pattern:'*.{fa,txt,log,yml}'

    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    tuple val(sampleName), path('*.scaffolds.fa'), path(read1), path(read2), emit: scaffolds
    path('*.log')    
    path "versions.yml"   

    script:
    """
    spades.py \\
        --threads $task.cpus \\
        --rnaviral \\
        -1 $read1 \\
        -2 $read2 \\
        -o ./

    mv spades.log ${sampleName}.spades.log
    if [ -f warnings.log ]; then
        mv warnings.log ${sampleName}.warnings.log
    fi    
    
    if [ -f scaffolds.fasta ]; then
        mv scaffolds.fasta ${sampleName}.scaffolds.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}