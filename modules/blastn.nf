process BLASTN {

    conda 'bioconda::blast=2.12.0'
    container 'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0'

    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.yml'

    input:
    tuple val(sampleName), path(scaffolds), path(read1), path(read2)
    path blast_db

    output:
    tuple val(sampleName), path("${sampleName}_blast.out"), path(scaffolds), path(read1), path(read2), path(blast_db), emit: blastn_out
    path "*.{log,sh,yml}"

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
