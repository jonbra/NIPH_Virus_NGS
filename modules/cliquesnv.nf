process CLIQUE_SNV {

  publishDir "${params.outdir}/2_clique_snv/", mode:'copy', pattern:'*.{fasta,json}'

  input:
  path samfile

  output:
  path ("*.fasta"), emit: CLIQUE_out
  path "*.json"

  script:
  """
  mkdir TEST
  #java -jar /usr/local/bin/clique-snv.jar -m snv-illumina -threads $task.cpu -outDir TEST -in ${samfile} -fdf extended
  cliquesnv -m snv-illumina -threads $task.cpu -outDir TEST -in ${samfile} -fdf extended
  """

}
