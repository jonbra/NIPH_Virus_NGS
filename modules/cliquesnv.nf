process CLIQUE_SNV {

  publishDir "${params.outdir}/2_clique_snv/", mode:'copy', pattern:'*.{fasta,json}'

  input:
  path samfile

  output:
  path ("TEST/*.fasta"), emit: CLIQUE_out
  path "TEST/*.json"

  script:
  """
  cliquesnv -m snv-illumina -threads $task.cpus -outDir TEST -in ${samfile} -fdf extended
  """

}
