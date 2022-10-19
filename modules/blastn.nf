process BLASTN {

    container 'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0'

    publishDir "${params.outdir}/5_blast/", mode:'copy', pattern:'*.{out,yml}'

    input:
    tuple val(sampleName), path(scaffolds)
    path blast_db

    output:
    tuple val(sampleName), path('*blast.out')     , optional:true, emit: blastn_out
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
2. Extract pairs of query and hit
*/