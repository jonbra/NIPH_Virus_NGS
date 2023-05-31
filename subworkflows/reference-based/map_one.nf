process MAP_ONE {
    
    container 'jonbra/viral_haplo:1.3'

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome
    path genome_index

    publishDir "${params.outdir}/ref-based", mode: 'copy', pattern:'*first_mapping*.{stats}'
    publishDir "${params.outdir}/ref-based", mode: 'copy', pattern:'*.{log,sh,txt,yml}'

    output:
    tuple val(sampleName), path ("${sampleName}.first_mapping.sorted.bam"), path ("${sampleName}.first_mapping.sorted.bam.bai"), optional: true, emit: sorted_out
    path "*.log", emit: BOWTIE2_log
    path "*.{stats,sh,txt}"

    script:
    """
    # First map against all references
    bowtie2 \
        --threads $task.cpus \
        --local \
        --very-sensitive-local \
        -x $genome \
        -1 ${read1} -2 ${read2} \
        2> ${sampleName}.bowtie2.first_mapping.log \
        -S ${sampleName}.first_mapping.sam
    samtools sort -@ $task.cpus -T ${sampleName} -O bam -o ${sampleName}.first_mapping.sorted.bam ${sampleName}.first_mapping.sam

    samtools index ${sampleName}.first_mapping.sorted.bam

    samtools stats ${sampleName}.first_mapping.sorted.bam > ${sampleName}.first_mapping.sorted.bam.stats

    SEQ=\$(grep 'reads mapped:' ${sampleName}.first_mapping.sorted.bam.stats| cut -f3)
    echo \${SEQ} > ${sampleName}.first_mapping_MAPPING_info.txt

    cp .command.sh ${sampleName}.bowtie2.first_mapping.sh

    cat <<-END_VERSIONS > bowtie2.first_mapping.versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
