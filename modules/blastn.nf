process BLASTN {

    container 'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0'

    publishDir "${params.outdir}/5_blast/", mode:'copy', pattern:'*.{yml}'

    input:
    tuple val(sampleName), path(scaffolds), path(read1), path(read2)
    path blast_db

    output:
    tuple val(sampleName), path('*blast.out'), path(scaffolds), path(read1), path(read2), emit: blastn_out
    path "versions.yml"   

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}

/*
Tasks:
1. Add header to blast out file - maybe not?
*/
