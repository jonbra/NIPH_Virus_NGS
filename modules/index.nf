process INDEX {

    container 'jonbra/viral_haplo:1.3'

    tag "$genome"
    publishDir "${params.outdir}/index", mode:'copy', pattern:'*.{bt2,log,sh,yml}'

    label 'small'

    input:
    path genome

    output:
    path "${genome}.*", emit: INDEX_out
    path "*.{log,sh}"

    script:
    """
    bowtie2-build --threads $task.cpus $genome $genome
    cp .command.log bowtie2_index.log
    cp .command.sh bowtie2_index.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
