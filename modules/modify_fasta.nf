process MODIFY_FASTA {

  tag "MODIFY_FASTA on $sampleName"

  publishDir "${params.outdir}/3_consensus", mode: 'copy', pattern:'*.fa'

  input:
  tuple val(sampleName), path(read1), path(read2)
  path consensus


  output:
  path "*.fa"

  script:
  """
  # Change the fasta header to include the sampleName
  bbrename.sh in=${consensus} out=${sampleName}.major_consensus.modified.fa prefix=${sampleName} addprefix=t
  """
}
