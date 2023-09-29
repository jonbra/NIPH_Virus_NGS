process BLASTN {

    conda 'bioconda::blast=2.12.0'
    container 'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0'

    publishDir "${params.outdir}/5_blast/", mode:'copy', pattern:'*.{out}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'

    input:
    tuple val(sampleName), path(scaffolds), path(read1), path(read2)
    path blast_db

    output:
    tuple val(sampleName), path("${sampleName}_blast.out"), path(scaffolds), path(read1), path(read2), path(blast_db), emit: blastn_out
    path("${sampleName}_blast.out")                                                                                  , emit: for_summarize
    path "*.{log,sh}"
    path "blastn_versions.yml", emit: versions

    script:

    """
    makeblastdb \\
        -in $blast_db \\
        -dbtype nucl

    blastn \\
        -db $blast_db \\
        -query $scaffolds \\
        -outfmt 6 \\
        -max_target_seqs 1 \\
        -out ${sampleName}_blast.out
    
    cp .command.log ${sampleName}.blastn_command.log
    cp .command.sh ${sampleName}.blastn_command.sh

    cat <<-END_VERSIONS > blastn_versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
