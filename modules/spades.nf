process SPADES {

    conda 'bioconda::spades=3.15.5'
    container 'quay.io/biocontainers/spades:3.15.4--h95f258a_0'

    // Force Nextflow to run a single sample at the time due to resource demands
    maxForks 1

    // Sometimes there can be zero reads after SUBSET_KRAKEN2, in which case Spades will crash.
    errorStrategy 'ignore'

    label 'large'

    publishDir "${params.outdir}/4_spades/", mode:'copy', pattern:'*.fa'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'

    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    tuple val(sampleName), path("${sampleName}.scaffolds.fa"), path(read1), path(read2), emit: scaffolds, optional: true
    path('*.log')     
    path "${sampleName}.spades_command.sh"  
    path "spades_versions.yml", emit: versions 

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