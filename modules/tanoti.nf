process TANOTI {
    tag "$sampleName"

    label 'large'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome

    publishDir "${params.outdir}/2_bam", mode: 'link', pattern:'*.{bam,bai}'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{stats,log,sh}'

    output:
    tuple val(sampleName), path ("${sampleName}.sorted.bam"), path ("${sampleName}.sorted.bam.bai"), emit: TANOTI_out
    path "*.log", emit: TANOTI_log
    path "${sampleName}_QC_PASS"
    path "*.{stats,sh}"

    script:
    """
    # Add mapper to path
    export PATH="~/Prosjekter/Tanoti/src:$PATH"

    #align vs. entire db
    tanoti \
        -r $genome \
        -i ${read1} ${read2} \
        -o ${sampleName}.sam \
        -p 1 -u 0 -m 85

    samtools view -bS ${sampleName}.sam | samtools sort bam -o ${sampleName}.sorted.bam

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
    tanoti \
        -r "\${major}" \
        -i ${read1} ${read2} \
        -o ${sampleName}.major.sam \
        -p 1 -u 0 -m 95

    samtools view -bS ${sampleName}.major.sam | samtools sort bam -o ${sampleName}.major.sorted.bam

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

    cp .command.sh ${sampleName}.tanoti.sh
    """
}
