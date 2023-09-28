process MAP_MINORITY_BOWTIE2 {

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
    tuple val(sampleName), path ("${sampleName}.minor.markdup.bam"), path ("${sampleName}.minor.markdup.bam.bai"), optional: true, emit: minority_out
    path "${sampleName}.minor.markdup.bam"                                                                       , optional: true, emit: GLUE
    path "*.log", optional: true, emit: BOWTIE2_log
    path "*.{stats,sh,txt}", optional: true

    script:
    """
    # Map reads against the reference with the second most mapped reads, regardless of the number of reads or coverage. 

    # Getting the reference with the second most mapped reads
    minor="\$(samtools idxstats ${sampleName}.first_mapping.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -2 | head -1 | cut -f1)" # Reference with most reads
    echo "\${minor}" >> ${sampleName}.minority_mapping_MAPPING_info.txt
    samtools faidx ${genome} "\${minor}" > "\${minor}".fa # Get the fasta sequence
    bowtie2-build "\${minor}".fa "\${minor}"
    bowtie2 \
        --threads $task.cpus \
        --local \
        --very-sensitive-local \
        -x "\${minor}" \
        -1 ${read1} -2 ${read2} \
        2> ${sampleName}.bowtie2.minority_mapping.log \
        | samtools sort -@ $task.cpus -T ${sampleName} -O bam - > ${sampleName}.minor.sorted.bam

    samtools index ${sampleName}.minor.sorted.bam

    samtools stats ${sampleName}.minor.sorted.bam > ${sampleName}.minor.sorted.bam.stats
    SEQ2=\$(grep 'reads mapped:' ${sampleName}.minor.sorted.bam | cut -f3)
    echo \${SEQ2} >> ${sampleName}.minority_mapping_MAPPING_info.txt

    # Then remove duplicates
    samtools sort -n ${sampleName}.minor.sorted.bam \
      | samtools fixmate -m - - \
      | samtools sort -O BAM \
      | samtools markdup --no-PG -r - ${sampleName}.minor.markdup.bam

    samtools index ${sampleName}.minor.markdup.bam

    samtools stats ${sampleName}.minor.markdup.bam > ${sampleName}.minor.markdup.bam.stats

    SEQ3=\$(grep 'reads mapped:' ${sampleName}.minor.markdup.bam.stats | cut -f3)
    echo \${SEQ3} >> ${sampleName}.minority_mapping_MAPPING_info.txt

    cp .command.sh ${sampleName}.bowtie2.minority_mapping.sh

    cat <<-END_VERSIONS > bowtie2.minority_mapping.versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}