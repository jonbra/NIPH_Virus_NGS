process DEDUP {

  tag "DEDUP on $sampleName"

  input:
  tuple val(sampleName), path(read1), path(read2)

  output:
  tuple val(sampleName), path ("${sampleName}_trimmed_dedup_R1.fastq"), path ("${sampleName}_trimmed_dedup_R2.fastq"), emit: DEDUP_out

  script:
  """
  # Change the fasta header to include the sampleName
  # Teh ac=f flag disables containment removal. I.e. only excat identicals will be removed.
  dedupe.sh threads=auto ac=f in1=${read1} in2=${read2} out=${sampleName}_trimmed_dedup.fastq
  reformat.sh in=${sampleName}_trimmed_dedup.fastq out1=${sampleName}_trimmed_dedup_R1.fastq out2=${sampleName}_trimmed_dedup_R2.fastq
  """
}
