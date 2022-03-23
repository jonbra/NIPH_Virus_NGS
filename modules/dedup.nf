process DEDUP {

  tag "DEDUP on $sampleName"

  input:
  tuple val(sampleName), path(read1), path(read2)

  output:
  tuple val(sampleName), path ("${sampleName}_trimmed_R1_dedup.fastq"), path ("${sampleName}_trimmed_R2_dedup.fastq"), emit: DEDUP_out

  script:
  """
  # Change the fasta header to include the sampleName
  dedupe.sh in=${read1} out=${sampleName}_trimmed_R1_dedup.fastq
  dedupe.sh in=${read2} out=${sampleName}_trimmed_R2_dedup.fastq
  """
}
