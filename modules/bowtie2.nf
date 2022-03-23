process BOWTIE {
    tag "$sampleName"
    errorStrategy 'ignore'

    label 'large'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome
    path genome_index

    publishDir "${params.outdir}/2_bam", mode: 'link', pattern:'*.{bam,bai}'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{stats,log,sh,txt}'

    output:
    tuple val(sampleName), path ("${sampleName}.major.sorted.bam"), path ("${sampleName}.major.sorted.bam.bai"), emit: BOWTIE2_out
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
        2> ${sampleName}.bowtie2.log \
        | samtools sort -@ $task.cpus -T ${sampleName} -O bam - > ${sampleName}.sorted.bam

    samtools index ${sampleName}.sorted.bam

    samtools stats ${sampleName}.sorted.bam > ${sampleName}.sorted.bam.stats

    SEQ=\$(grep 'reads mapped:' ${sampleName}.sorted.bam.stats | cut -f3)
    echo \${SEQ} > ${sampleName}_MAPPING_info.txt

    # Then re-map against the reference with the most mapped reads
    major="\$(samtools idxstats ${sampleName}.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)" # Reference with most reads
    echo "\${major}" >> ${sampleName}_MAPPING_info.txt
    samtools faidx ${genome} "\${major}" > "\${major}".fa # Get the fasta sequence
    bowtie2-build "\${major}".fa "\${major}"
    bowtie2 \
        --threads $task.cpus \
        --local \
        --very-sensitive-local \
        -x "\${major}" \
        -1 ${read1} -2 ${read2} \
        2> ${sampleName}.bowtie2_major.log \
        | samtools sort -@ $task.cpus -T ${sampleName} -O bam - > ${sampleName}.major.sorted.bam

    samtools index ${sampleName}.major.sorted.bam

    samtools stats ${sampleName}.major.sorted.bam > ${sampleName}.major.sorted.bam.stats
    SEQ2=\$(grep 'reads mapped:' ${sampleName}.major.sorted.bam.stats | cut -f3)
    echo \${SEQ2} >> ${sampleName}_MAPPING_info.txt


    cp .command.sh ${sampleName}.bowtie2.sh
    """
}
/*
if [[ \$SEQ -gt 999 ]]
then
    echo 'More than 1000 reads in first mapping' > ${sampleName}_MAPPING_info.txt

    # Then re-map against the reference with the most mapped reads
    major="\$(samtools idxstats ${sampleName}.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)" # Reference with most reads
    echo "\${major}" >> ${sampleName}_MAPPING_info.txt
    samtools faidx ${genome} "\${major}" > "\${major}".fa # Get the fasta sequence
    bowtie2-build "\${major}".fa "\${major}"
    bowtie2 \
        --threads $task.cpus \
        --local \
        --very-sensitive-local \
        -x "\${major}" \
        -1 ${read1} -2 ${read2} \
        2> ${sampleName}.bowtie2_major.log \
        | samtools sort -@ $task.cpus -T ${sampleName} -O bam - > ${sampleName}.major.sorted.bam

    samtools index ${sampleName}.major.sorted.bam

    samtools stats ${sampleName}.major.sorted.bam > ${sampleName}.major.sorted.bam.stats

    # Remove duplicates
    samtools rmdup ${sampleName}.major.sorted.bam ${sampleName}.major.sorted.rmdup.bam
    samtools index ${sampleName}.major.sorted.rmdup.bam

    samtools stats ${sampleName}.major.sorted.bam > ${sampleName}.major.sorted.bam.stats

else
    echo 'Less than 1000 reads in first mapping' > ${sampleName}_MAPPING_info.txt
fi
*/
