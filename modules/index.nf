process INDEX {

    conda 'bioconda::bowtie2=2.4.4'
    container 'quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0'

    publishDir "${params.outdir}/index", mode:'copy', pattern:'*.{bt2,log,sh}'

    label 'small'

    input:
    path genome

    output:
    path "${genome}.*", emit: INDEX_out
    path "*.{log,sh}"
    path "index_versions.yml", emit: versions

    script:
    """
    bowtie2-build --threads $task.cpus $genome $genome
    cp .command.log bowtie2_index.log
    cp .command.sh bowtie2_index.sh

    cat <<-END_VERSIONS > index_versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
