process ABACAS {

    container 'quay.io/biocontainers/abacas:1.3.1--pl526_0'

    publishDir "${params.outdir}/6_abacas/", mode:'copy', pattern:'*'

    input:
    path reference 
    path scaffolds
    tuple val(sampleName), path('*.tsv')

    output:
    tuple val(sampleName), path('*.abacas*'), emit: abacas_results
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    for i in \$(ls *scaffolds.fa)
    do
    # Extract the first characters before underscore
    GENO=\$(echo \$i | cut -d'_' -f 1)
    # Pull out the reference
    REF=\$(ls \$GENO*ref.fa)
    # Pull out the coresponding scaffolds
    SCAF=\$(ls \$GENO*scaffolds.fa)

    # Run ABACAS
    abacas.pl \\
        -r \$REF \\
        -q \$SCAF \\
        -p nucmer \\
        -o ${sampleName}_\$GENO.abacas
    mv nucmer.delta ${sampleName}.abacas.nucmer.delta
    mv nucmer.filtered.delta ${sampleName}.abacas.nucmer.filtered.delta
    mv nucmer.tiling ${sampleName}.abacas.nucmer.tiling
    mv unused_contigs.out ${sampleName}.abacas.unused.contigs.out
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abacas: \$(echo \$(abacas.pl -v 2>&1) | sed 's/^.*ABACAS.//; s/ .*\$//')
    END_VERSIONS
    """
}