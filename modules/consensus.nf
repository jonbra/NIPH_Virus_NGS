process CONSENSUS {
  tag "CONSENSUS on $sampleName"
  errorStrategy 'ignore'

  publishDir "${params.outdir}/3_consensus", mode: 'copy', pattern:'*.{fa,csi,bed,gz}'
  publishDir "${params.outdir}/3_consensus/log", mode: 'link', pattern:'*.{stats,log,sh}'

  input:
  tuple val(sampleName), path(read1), path(read2)
  tuple val(sampleName), path ("${sampleName}.major.sorted.bam"), path ("${sampleName}.major.sorted.bam.bai")
  path genome

  output:
  path "*.major_consensus.fa", emit: CONSENSUS_fa
  path "*.{gz,csi,bed}"


  """
  # Make the major genotype again
  mapped_reads="\$(samtools idxstats ${sampleName}.major.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f2)"
  if [ "\${mapped_reads}" -gt 499 ]
  then
    major="\$(samtools idxstats ${sampleName}.major.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)" # Reference with most reads
    samtools faidx ${genome} "\${major}" > "\${major}".fa # Get the fasta sequence

    # Order and fix the bam file
    samtools sort -n ${sampleName}.major.sorted.bam > ${sampleName}.major.sorted.byQuery.bam
    samtools fixmate -m ${sampleName}.major.sorted.byQuery.bam ${sampleName}.major.sorted.fix.bam
    samtools sort ${sampleName}.major.sorted.fix.bam > ${sampleName}.major.sorted.fix.sorted.bam

    # Make the consensus using the best reference
    samtools markdup -r ${sampleName}.major.sorted.fix.sorted.bam ${sampleName}.major.sorted.fix.sorted.marked.bam
    bcftools mpileup -f "\${major}".fa ${sampleName}.major.sorted.fix.sorted.marked.bam | bcftools call -mv -Ob -o ${sampleName}.calls.vcf.gz
    bcftools index ${sampleName}.calls.vcf.gz
    samtools index ${sampleName}.major.sorted.fix.sorted.marked.bam

    # Get the coverage for each position
    bedtools genomecov -bga -ibam ${sampleName}.major.sorted.fix.sorted.marked.bam| awk '\$4 < 6' > regionswithlessthan6coverage.bed
    bcftools consensus -m regionswithlessthan6coverage.bed -f "\${major}".fa calls.vcf.gz -o ${sampleName}.major_consensus.fa
  else
  echo "Fewer than 500 mapped reads" >> ${sampleName}.consensus.log
  fi
  """
}
