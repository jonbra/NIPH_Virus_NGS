process CONSENSUS {
  tag "CONSENSUS on $sampleName"

  publishDir "${params.outdir}/3_consensus", mode: 'copy', pattern:'*.{fa,bam,gz,csi,bed,bai}'
  publishDir "${params.outdir}/3_consensus/log", mode: 'link', pattern:'*.{stats,log,sh}'

  input:
  tuple val(sampleName), path ("${sampleName}.sorted.bam"), path ("${sampleName}.sorted.bam.bai")
  path ref_file

  output:
  tuple val(sampleName), path ("${sampleName}.sorted.fix.sorted.marked.bam"), emit: CONSENSUS_out
  path "*.{fa,gz,csi,bed,bai}"

  script:
  """
  # First get the name of the best reference genotype
  major="\$(samtools idxstats ${sampleName}.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)"

  # Then get the fasta sequence of the best reference
  samtools faidx ${ref_file} "\${major}" > "\${major}".fa

  # Order and fix the bam file
  samtools sort -n ${sampleName}.sorted.bam > ${sampleName}.sorted.byQuery.bam
  samtools fixmate -m ${sampleName}.sorted.byQuery.bam ${sampleName}.sorted.fix.bam
  samtools sort ${sampleName}.sorted.fix.bam > ${sampleName}.sorted.fix.sorted.bam

  # Make the consensus using the best reference
  samtools markdup -r ${sampleName}.sorted.fix.sorted.bam ${sampleName}.sorted.fix.sorted.marked.bam
  bcftools mpileup -f "\${major}".fa ${sampleName}.sorted.fix.sorted.marked.bam | bcftools call -mv -Ob -o calls.vcf.gz
  bcftools index calls.vcf.gz
  samtools index ${sampleName}.sorted.fix.sorted.marked.bam

  # Get the coverage for each position
  bedtools genomecov -bga -ibam ${sampleName}.sorted.fix.sorted.marked.bam | awk '\$4 < 6' > regionswithlessthan6coverage.bed
  bcftools consensus -m regionswithlessthan6coverage.bed -f "\${major}".fa calls.vcf.gz -o ${sampleName}_consensus.fa
  """
}
/* GJENSTÅR
Denne kan erstattes med en sed kommando:
seqkit replace -p "(.+)" -r ${bestF3%%_*} cons.fa > ${bestF3%%_*}_consensus.fa #endrer navn fra referanse-navn til prøvenavn inne i fasta-fil
*/
