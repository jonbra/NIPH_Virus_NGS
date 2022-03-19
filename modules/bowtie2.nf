process BOWTIE {
    tag "$sampleName"
    errorStrategy 'ignore'

    label 'large'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome
    path genome_index

    publishDir "${params.outdir}/2_bam", mode: 'link', pattern:'*.{bam,bai}'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{stats,log,sh}'

    output:
    tuple val(sampleName), path ("${sampleName}.major.sorted.bam"), path ("${sampleName}.major.sorted.bam.bai"), emit: BOWTIE2_out
    path "*.log", emit: BOWTIE2_log
    path "${sampleName}_QC_*"
    path "*.{stats,sh}"

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

      SEQ=\$(grep 'raw total sequences against all references:' ${sampleName}.sorted.bam.stats | cut -f3)
      if [[ \$SEQ -gt 1000 ]]
      then
          echo 'PASS' > ${sampleName}_QC_PASS
      else
          echo 'FAIL' > ${sampleName}_QC_FAIL
      fi

    # Then re-map against the reference with the most mapped reads
    major="\$(samtools idxstats ${sampleName}.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)" # Reference with most reads
    samtools faidx ${genome} "\${major}" > "\${major}".fa # Get the fasta sequence
    bowtie2-build "\${major}".fa "\${major}"
    bowtie2 \
        --threads $task.cpus \
        --local \
        --very-sensitive-local \
        -x "\${major}" \
        -1 ${read1} -2 ${read2} \
        -S ${sampleName}.major.bam \
        2> ${sampleName}.bowtie2_major.log \
        | samtools sort -@ $task.cpus -T ${sampleName} -O bam - > ${sampleName}.major.sorted.bam

    if test -f "${sampleName}.major.sorted.bam";
    then
      samtools index ${sampleName}.major.sorted.bam

      samtools stats ${sampleName}.major.sorted.bam > ${sampleName}.major.sorted.bam.stats

      SEQ=\$(grep 'raw total sequences against best reference:' ${sampleName}.major.sorted.bam.stats | cut -f3)
    else
      echo "${sampleName}.major.sorted.bam" does not exist
    fi

    if [[ \$SEQ -gt 1000 ]]
    then
        echo 'PASS' >> ${sampleName}_QC_PASS
    else
        echo 'FAIL' >> ${sampleName}_QC_FAIL
    fi

    cp .command.sh ${sampleName}.bowtie2.sh

    """
}
