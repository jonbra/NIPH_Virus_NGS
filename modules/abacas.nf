process ABACAS {

    container 'quay.io/biocontainers/abacas:1.3.1--pl526_0'
    
    errorStrategy 'ignore'

    publishDir "${params.outdir}/6_abacas/", mode:'copy', pattern:'*'

    input:
    tuple val(sampleName), path(references), path(scaffolds)

    output:
    tuple val(sampleName), path('*.abacas*'), emit: abacas_results
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    # TODO: Få inn samplename også i fasta-headerne
    
    for i in \$(ls *scaffolds.fa)
    do
        # Extract the genotype (first cut on underscore, then on dot)
        GENO=\$(echo \$i | cut -d '_' -f 1 | cut -d '.' -f 2)
        # Extract the first field (samplename) when cut on dot
        # SAMPLE=\$(echo \$i | cut -d',' -f 1)
        # Pull out the reference
        REF=\$(ls *\${GENO}_ref.fa)
        # Pull out the corresponding scaffolds
        SCAF=\$(ls *\${GENO}_scaffolds.fa)

        # Run ABACAS
        abacas.pl \\
            -r \$REF \\
            -q \$SCAF \\
            -p nucmer \\
            -m \\
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