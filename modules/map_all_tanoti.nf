process MAP_ALL_TANOTI {
    
    container 'jonbra/viral_haplo:1.3'

    label 'medium'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path references

    publishDir "${params.outdir}/6_map", mode: 'copy', pattern:'*first_mapping*.{stats}'
    publishDir "${params.outdir}/6_map", mode: 'copy', pattern:'*.{log,sh,txt,yml}'

    output:
    tuple val(sampleName), path ("${sampleName}.first_mapping.sorted.bam"), path ("${sampleName}.first_mapping.sorted.bam.bai"), optional: true, emit: sorted_out
    path "*.log"                                                                                                               , optional: true, emit: TANOTI_log
    path "*.{stats,sh,txt}"                                                                                                    , optional: true

    script:
    """
    # First map against all references
    tanoti \
        -r $references \
        -i ${read1} ${read2} \
        -o ${sampleName}.first_mapping.sam \
        -p 1 -u 1 -P 8 -m ${params.tanoti_stringency_1}

    samtools sort -@ $task.cpus -T ${sampleName} -O bam -o ${sampleName}.first_mapping.sorted.bam ${sampleName}.first_mapping.sam
    samtools index ${sampleName}.first_mapping.sorted.bam
    samtools stats ${sampleName}.first_mapping.sorted.bam > ${sampleName}.first_mapping.sorted.bam.stats

    SEQ=\$(grep 'reads mapped:' ${sampleName}.first_mapping.sorted.bam.stats| cut -f3)
    echo \${SEQ} > ${sampleName}.first_mapping_MAPPING_info.txt

    cp .command.sh ${sampleName}.tanoti.first_mapping.sh

    cat <<-END_VERSIONS > tanoti.first_mapping.versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

 