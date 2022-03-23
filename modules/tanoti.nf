process TANOTI {
    tag "$sampleName"
    errorStrategy 'ignore'

    label 'large'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome

    publishDir "${params.outdir}/2_bam", mode: 'copy', pattern:'*.{bam,bai}'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{stats,sh}'

    output:
    tuple val(sampleName), path ("${sampleName}.major.sorted.fix.sorted.marked.bam"), path ("${sampleName}.major.sorted.fix.sorted.marked.bam.bai"), emit: TANOTI_out
    path "${sampleName}_QC_*"
    path "*.{stats,sh}"

    script:
    """
    #align vs. entire db
    tanoti \
        -r $genome \
        -i ${read1} ${read2} \
        -o ${sampleName}.sam \
        -p 1 -u 0 -m 85

    samtools view -bS ${sampleName}.sam | samtools sort -o ${sampleName}.sorted.bam

    samtools index ${sampleName}.sorted.bam

    samtools stats ${sampleName}.sorted.bam > ${sampleName}.sorted.bam.stats

    SEQ=\$(grep 'raw total sequences' ${sampleName}.sorted.bam.stats | cut -f3)
    if [[ \$SEQ -gt 999 ]]
    then
        echo 'PASS' > ${sampleName}_QC_MAPPING_PASS
        # Re-map against the reference with the most mapped reads
        major="\$(samtools idxstats ${sampleName}.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)" # Reference with most reads
        samtools faidx ${genome} "\${major}" > "\${major}".fa # Get the fasta sequence
        tanoti \
            -r "\${major}".fa \
            -i ${read1} ${read2} \
            -o ${sampleName}.major.sam \
            -p 1 -u 0 -m 95
        samtools view -bS ${sampleName}.major.sam | samtools sort -o ${sampleName}.major.sorted.bam

        samtools index ${sampleName}.major.sorted.bam

        samtools stats ${sampleName}.major.sorted.bam > ${sampleName}.major.sorted.bam.stats

        # Order and fix the bam file
        samtools sort -n ${sampleName}.major.sorted.bam > ${sampleName}.major.sorted.byQuery.bam
        samtools fixmate -m ${sampleName}.major.sorted.byQuery.bam ${sampleName}.major.sorted.fix.bam
        samtools sort ${sampleName}.major.sorted.fix.bam > ${sampleName}.major.sorted.fix.sorted.bam

        # Remove duplicates
        samtools markdup -r ${sampleName}.major.sorted.fix.sorted.bam ${sampleName}.major.sorted.fix.sorted.marked.bam
        samtools index ${sampleName}.major.sorted.fix.sorted.marked.bam

        samtools stats ${sampleName}.major.sorted.fix.sorted.marked.bam > ${sampleName}.major.sorted.fix.sorted.marked.bam.stats  
    else
        echo 'FAIL' > ${sampleName}_QC_MAPPING_FAIL
    fi

    cp .command.sh ${sampleName}.tanoti.sh
    """
}
