process MAP_MAJORITY {

    container 'jonbra/viral_haplo:1.3'

    label 'small'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome
    path genome_index
    tuple val(sampleName), path ("${sampleName}.first_mapping.sorted.bam"), path ("${sampleName}.first_mapping.sorted.bam.bai")

    publishDir "${params.outdir}/ref-based", mode: 'copy', pattern:'*major*.{bam,bai}'
    publishDir "${params.outdir}/ref-based", mode: 'copy', pattern:'*.{stats,log,sh,txt,yml}'

    output:
    //tuple val(sampleName), path ("${sampleName}.markdup.bam"), path ("${sampleName}.markdup.bam.bai"), optional: true, emit: markdup_out
    tuple val(sampleName), path ("${sampleName}.major.markdup.bam"), path ("${sampleName}.major.markdup.bam.bai"), optional: true, emit: majority_out
    path "${sampleName}.major.markdup.bam"                                                                       , emit: GLUE
    path "*.log", emit: BOWTIE2_log
    path "*.{stats,sh,txt}"

    script:
    """
    # Re-map against the reference with the most mapped reads
    major="\$(samtools idxstats ${sampleName}.first_mapping.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)" # Reference with most reads
    echo "\${major}" >> ${sampleName}.majority_mapping_MAPPING_info.txt
    samtools faidx ${genome} "\${major}" > "\${major}".fa # Get the fasta sequence
    bowtie2-build "\${major}".fa "\${major}"
    bowtie2 \
        --threads $task.cpus \
        --local \
        --very-sensitive-local \
        -x "\${major}" \
        -1 ${read1} -2 ${read2} \
        2> ${sampleName}.bowtie2.majority_mapping.log \
        | samtools sort -@ $task.cpus -T ${sampleName} -O bam - > ${sampleName}.major.sorted.bam

    samtools index ${sampleName}.major.sorted.bam

    samtools stats ${sampleName}.major.sorted.bam > ${sampleName}.major.sorted.bam.stats
    SEQ2=\$(grep 'reads mapped:' ${sampleName}.major.sorted.bam.stats | cut -f3)
    echo \${SEQ2} >> ${sampleName}.majority_mapping_MAPPING_info.txt

    # Then remove duplicates
    samtools sort -n ${sampleName}.major.sorted.bam \
      | samtools fixmate -m - - \
      | samtools sort -O BAM \
      | samtools markdup --no-PG -r - ${sampleName}.major.markdup.bam

    samtools index ${sampleName}.major.markdup.bam

    samtools stats ${sampleName}.major.markdup.bam > ${sampleName}.major.markdup.bam.stats

    SEQ3=\$(grep 'reads mapped:' ${sampleName}.major.markdup.bam.stats | cut -f3)
    echo \${SEQ3} >> ${sampleName}.majority_mapping_MAPPING_info.txt

    cp .command.sh ${sampleName}.bowtie2.majority_mapping.sh

    cat <<-END_VERSIONS > bowtie2.majority_mapping.versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
