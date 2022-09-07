process ABACAS {

    publishDir "${params.outdir}/abacas/", mode:'copy', pattern:'*'

    input:
    tuple val(sampleName), path(scaffolds)
    path ref_file

    output:
    tuple val(sampleName), path('*.abacas*'), emit: abacas_results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    abacas.pl \\
        -r $ref_file \\
        -q $scaffolds \\
        -p nucmer \\
        -o ${sampleName}.abacas
    mv nucmer.delta ${sampleName}.abacas.nucmer.delta
    mv nucmer.filtered.delta ${sampleName}.abacas.nucmer.filtered.delta
    mv nucmer.tiling ${sampleName}.abacas.nucmer.tiling
    mv unused_contigs.out ${sampleName}.abacas.unused.contigs.out
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abacas: \$(echo \$(abacas.pl -v 2>&1) | sed 's/^.*ABACAS.//; s/ .*\$//')
    END_VERSIONS
    """
}

/*
Thoughts:
How to separate out the different references and run abacas for each of them?
*/