process CONSENSUS {
  tag "CONSENSUS on $sampleName"

  publishDir "${params.outdir}/3_consensus", mode: 'copy', pattern:'*.{fa,bam,gz,csi,bed,bai}'
  publishDir "${params.outdir}/3_consensus/log", mode: 'link', pattern:'*.{stats,log,sh}'

  input:
  tuple val(sampleName), path(read1), path(read2)
  tuple val(sampleName), path ("${sampleName}.sorted.bam"), path ("${sampleName}.sorted.bam.bai")
  path ref_file

  output:
  tuple val(sampleName), path ("${sampleName}.major.sorted.fix.sorted.marked.bam")
  path "*.major_consensus.fa", emit: CONSENSUS_fa
  path "*.{gz,csi,bed,bai}"


  """
  # First get the name of the best reference genotype
  major="\$(samtools idxstats ${sampleName}.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)"

  # Then get the fasta sequence of the best reference
  samtools faidx ${ref_file} "\${major}" > "\${major}".fa

  # Need to remap here
  bowtie2-build "\${major}".fa "\${major}"
  bowtie2 \
      --threads $task.cpus \
      --local \
      --very-sensitive-local \
      -x "\${major}" \
      -1 ${read1} -2 ${read2} \
      -S ${sampleName}.major.bam \
      2> ${sampleName}.bowtie2_major.log

  # Sort and index
  samtools view -bS ${sampleName}.major.bam | samtools sort -o ${sampleName}.major.sorted.bam
  samtools index ${sampleName}.major.sorted.bam

  # Order and fix the bam file
  samtools sort -n ${sampleName}.major.sorted.bam > ${sampleName}.major.sorted.byQuery.bam
  samtools fixmate -m ${sampleName}.major.sorted.byQuery.bam ${sampleName}.major.sorted.fix.bam
  samtools sort ${sampleName}.major.sorted.fix.bam > ${sampleName}.major.sorted.fix.sorted.bam

  # Make the consensus using the best reference
  samtools markdup -r ${sampleName}.major.sorted.fix.sorted.bam ${sampleName}.major.sorted.fix.sorted.marked.bam
  bcftools mpileup -f "\${major}".fa ${sampleName}.major.sorted.fix.sorted.marked.bam | bcftools call -mv -Ob -o calls.vcf.gz
  bcftools index calls.vcf.gz
  samtools index ${sampleName}.major.sorted.fix.sorted.marked.bam

  # Get the coverage for each position
  bedtools genomecov -bga -ibam ${sampleName}.major.sorted.fix.sorted.marked.bam| awk '\$4 < 6' > regionswithlessthan6coverage.bed
  bcftools consensus -m regionswithlessthan6coverage.bed -f "\${major}".fa calls.vcf.gz -o ${sampleName}.major_consensus.fa
  """
}
