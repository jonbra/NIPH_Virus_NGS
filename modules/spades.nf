process SPADES {

    // Sometimes there can be zero reads after SUBSET_KRAKEN2, in which case Spades will crash.
    // Perhaps better to use an if statement to only execute the process in case there are reads in the input file
    errorStrategy 'ignore'

    publishDir "${params.outdir}/spades/", mode:'copy', pattern:'*.{fa,txt,log,yml}'

    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    tuple val(sampleName), path('*.scaffolds.fa')     , optional:true, emit: scaffolds
    tuple val(sampleName), path('*.log')                , emit: log
    path "versions.yml"   
    
    when:
    task.ext.when == null || task.ext.when

    script:
    """
    spades.py \\
        --threads $task.cpus \\
        --rnaviral \\
        -1 $read1 \\
        -2 $read2 \\
        -o ./

    mv spades.log ${sampleName}.spades.log
    
    if [ -f scaffolds.fasta ]; then
        mv scaffolds.fasta ${sampleName}.scaffolds.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}