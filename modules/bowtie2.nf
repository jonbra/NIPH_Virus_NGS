process BOWTIE2 {

    container 'jonbra/viral_haplo:1.2'

    tag "$sampleName"
    //errorStrategy 'ignore'

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome
    path genome_index

    publishDir "${params.outdir}/2_bam", mode: 'copy', pattern:'*.{bam,bai}'
    publishDir "${params.outdir}/2_bam", mode: 'copy', pattern:'*.{stats,log,sh,txt,yml}'

    output:
    //tuple val(sampleName), path ("${sampleName}.markdup.bam"), path ("${sampleName}.markdup.bam.bai"), optional: true, emit: markdup_out
    tuple val(sampleName), path ("${sampleName}.sorted.bam"), path ("${sampleName}.sorted.bam.bai"), optional: true, emit: sorted_out
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
        -S ${sampleName}.sam
    samtools sort -@ $task.cpus -T ${sampleName} -O bam -o ${sampleName}.sorted.bam ${sampleName}.sam

    samtools index ${sampleName}.sorted.bam

    samtools stats ${sampleName}.sorted.bam > ${sampleName}.sorted.bam.stats

    SEQ=\$(grep 'reads mapped:' ${sampleName}.sorted.bam.stats | cut -f3)
    echo \${SEQ} > ${sampleName}_MAPPING_info.txt

        cp .command.sh ${sampleName}.bowtie2.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

/*
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

    # Then remove duplicates
    samtools sort -n ${sampleName}.major.sorted.bam \
      | samtools fixmate -m - - \
      | samtools sort -O BAM \
      | samtools markdup --no-PG -r - ${sampleName}.markdup.bam

    samtools index ${sampleName}.markdup.bam

    samtools stats ${sampleName}.markdup.bam > ${sampleName}.markdup.bam.stats

    SEQ3=\$(grep 'reads mapped:' ${sampleName}.markdup.bam.stats | cut -f3)
    echo \${SEQ3} >> ${sampleName}_MAPPING_info.txt
*/
