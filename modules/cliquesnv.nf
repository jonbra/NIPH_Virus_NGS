process CLIQUE_SNV {
  publishDir "${params.outdir}/clique_snv/", mode:'copy', pattern:'*.{zip}'

  JEG KAN KLONE REPO OG KJØRE MAKE INSTALL

  input:
  tuple val(sampleName), path(read1), path(read2)
  val(source)

  output:
  tuple val(sampleName), path ("*.zip"), emit: FASTQC_out
  path "*.{log,sh}"

  script:
  """
  fastqc -t $task.cpus ${read1} ${read2}
  cp .command.log ${sampleName}.${source}.fastqc.log
  cp .command.sh ${sampleName}.${source}.fastqc.sh
  """

}
