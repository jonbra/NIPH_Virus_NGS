process SPADES {

    container 'quay.io/biocontainers/spades:3.15.4--h95f258a_0'

    maxForks 1

    // Sometimes there can be zero reads after SUBSET_KRAKEN2, in which case Spades will crash.
    errorStrategy 'ignore'

    label 'spades'

    publishDir "${params.outdir}/4_spades/", mode:'copy', pattern:'*.fa'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.yml'

    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    tuple val(sampleName), path("${sampleName}.scaffolds.fa"), path(read1), path(read2), emit: scaffolds
    path('*.log')    
    path "spades_versions.yml"   
    path "${sampleName}.spades_command.sh"   

    script:
    """
    spades.py \\
        --threads $task.cpus \\
        --$params.spades_mode \\
        -1 $read1 \\
        -2 $read2 \\
        -o ./

    mv spades.log ${sampleName}.spades.log
    if [ -f warnings.log ]; then
        mv warnings.log ${sampleName}.spades_warnings.log
    fi    
    
    if [ -f scaffolds.fasta ]; then
        mv scaffolds.fasta ${sampleName}.scaffolds.fa
    fi

    cp .command.log ${sampleName}.spades_command.log
    cp .command.sh ${sampleName}.spades_command.sh
    
    cat <<-END_VERSIONS > spades_versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}