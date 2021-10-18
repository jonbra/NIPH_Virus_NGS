process INDEX {
    tag "$genome"
    publishDir "${params.outdir}/index", mode:'copy', pattern:'*.bt2'

    label 'small'

    input:
    path genome

    output:
    path "${genome}.*", emit: INDEX_out
    path "*.{log,sh}"

    script:
    """
    bowtie2-build --seed 1 $genome $genome
    cp .command.log bowtie2_index.log
    cp .command.sh bowtie2_index.sh
    """
}
