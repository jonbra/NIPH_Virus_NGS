process REPAIR {

  container 'quay.io/biocontainers/bbmap:38.86--h1296035_0'

  errorStrategy 'terminate'

  publishDir "${params.outdir}/3_kraken2/", mode:'copy', pattern:'*.fastq'
  publishDir "${params.outdir}/3_kraken2/log/", mode:'copy', pattern:'*.{sh,out}'

  input:
  tuple val(sampleName), path(read1), path(read2)

  output:
  tuple val(sampleName), path ("*R1_repair.fastq"), path ("*R2_repair.fastq"), emit: paired_reads
  path "*.{sh,out}"

  script:
  """
  repair.sh -Xmx8g in=${read1} in2=${read2} out1=${sampleName}_R1_repair.fastq  out2=${sampleName}_R2_repair.fastq
  cp .command.sh ${sampleName}.repair.sh
  cp .command.out ${sampleName}.repair.out
  """
}