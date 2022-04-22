process CONSENSUS {
  tag "CONSENSUS on $sampleName"


  publishDir "${params.outdir}/3_consensus", mode: 'copy', pattern:'*.{fa,csi,bed,gz}'
  publishDir "${params.outdir}/3_consensus/log", mode: 'link', pattern:'*.{stats,txt,sh}'

  input:
  tuple val(sampleName), path ("${sampleName}.markdup.bam"), path ("${sampleName}.markdup.bam.bai")
  path genome

  output:
  path "*.major_consensus.fa", emit: CONSENSUS_fa
  path "*.{gz,csi,bed}"
  path "*.txt"


  """
  samtools stats ${sampleName}.markdup.bam > ${sampleName}.markdup.bam.stats
  SEQ=\$(grep 'reads mapped:' ${sampleName}.markdup.bam.stats | cut -f3)

  if [[ \$SEQ -gt 499 ]]
  then
    echo 'More than 500 reads (duplicates removed) mapped to best reference'> ${sampleName}_CONSENSUS_info.txt
    # Make the consensus using the reference with most mapped reads
    major="\$(samtools idxstats ${sampleName}.markdup.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)" # Reference with most reads
    samtools faidx ${genome} "\${major}" > "\${major}".fa # Get the fasta sequence
    bcftools mpileup -f "\${major}".fa ${sampleName}.markdup.bam | bcftools call -mv -Ob -o ${sampleName}.calls.vcf.gz
    bcftools index ${sampleName}.calls.vcf.gz
    samtools index ${sampleName}.markdup.bam

    # Get the coverage for each position
    bedtools genomecov -bga -ibam ${sampleName}.markdup.bam | awk '\$4 < 6' > regionswithlessthan6coverage.bed
    bcftools consensus -m regionswithlessthan6coverage.bed -f "\${major}".fa ${sampleName}.calls.vcf.gz -o ${sampleName}.major_consensus.fa

  else
    echo 'Less than 500 reads (duplicates removed) mapped to best reference'> ${sampleName}_CONSENSUS_info.txt
  fi
  """
}
