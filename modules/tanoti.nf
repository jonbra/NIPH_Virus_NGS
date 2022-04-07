process TANOTI {
    tag "$sampleName"

    label 'large'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome

    publishDir "${params.outdir}/2_bam", mode: 'copy', pattern:'*.{bam,bai}'
    publishDir "${params.outdir}/2_bam/log", mode: 'copy', pattern:'*.{stats,sh,txt}'

    output:
    tuple val(sampleName), path ("${sampleName}.major.sorted.bam"), path ("${sampleName}.major.sorted.bam.bai"), emit: TANOTI_out
    path "*.{stats,sh,txt}"

    script:
    """
    # First map against all references
    tanoti \
        -r $genome \
        -i ${read1} ${read2} \
        -o ${sampleName}.sam \
        -p 1 -u 0 -m 85

    samtools view -bS ${sampleName}.sam | samtools sort -o ${sampleName}.sorted.bam

    samtools index ${sampleName}.sorted.bam

    samtools stats ${sampleName}.sorted.bam > ${sampleName}.sorted.bam.stats

    SEQ=\$(grep 'reads mapped:' ${sampleName}.sorted.bam.stats | cut -f3)
    echo \${SEQ} > ${sampleName}_MAPPING_info.txt

    # Then re-map against the reference with the most mapped reads
    major="\$(samtools idxstats ${sampleName}.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)" # Reference with most reads
    echo "\${major}" >> ${sampleName}_MAPPING_info.txt
    samtools faidx ${genome} "\${major}" > "\${major}".fa # Get the fasta sequence

    tanoti \
      -r "\${major}".fa \
      -i ${read1} ${read2} \
      -o ${sampleName}.major.sam \
      -p 1 -u 0 -m 95

    samtools view -bS ${sampleName}.major.sam | samtools sort -o ${sampleName}.major.sorted.bam

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

    cp .command.sh ${sampleName}.tanoti.sh
    """
}
