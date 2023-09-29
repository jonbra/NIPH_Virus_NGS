process HBV_RT_BLAST {

    container 'quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0'

    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.yml'

    input:
    path 'scaffolds/'
    path reference

    output:
    path "*.json", emit: rt_blast
    path "*.{log,sh,yml}"

    script:

    """
    cd scaffolds/
    for i in \$(ls *scaffolds.fa)
    do

    #TODO: Determine the genotype of the query sequence 
    # and then choose subject sequence accordingly. 
    # Can be done with an if statement inside the loop?

    blastx \\
        -query \${i} \\
        -subject ../$reference \\
        -outfmt 15 \\
        -qcov_hsp_perc 20 > ../\${i%scaffolds.fa}_rt_blastx.json

    done

    cd ../
    cp .command.log blastx_command.log
    cp .command.sh blastx_command.sh

    cat <<-END_VERSIONS > blastx_versions.yml
    "${task.process}":
        blast: \$(blastx -version 2>&1 | sed 's/^.*blastx: //; s/ .*\$//')
    END_VERSIONS
    """
}
