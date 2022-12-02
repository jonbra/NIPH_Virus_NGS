process RVA_GENO {

    container 'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0'

    publishDir "${params.outdir}/7_genotyping/", mode:'copy', pattern:'*.{out,yml}'

    input:
    tuple val(sampleName), path(scaffolds)
    path genotypes

    output:
    tuple val(sampleName), path('*blast.out')     , optional:true, emit: genotyping_out
    path "versions.yml"   

    script:

    """
    makeblastdb \\
        -in $genotypes \\
        -dbtype nucl

    blastn \\
        -db $genotypes \\
        -query $scaffolds \\
        -outfmt 6 \\
        -max_target_seqs 5 \\
        -out ${sampleName}_blast.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}