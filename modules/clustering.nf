process CLUSTER {

  publishDir "${params.outdir}/2_clustering/", mode:'copy', pattern:'*.{treefile,aligned.fasta,log}'

  input:
  path 'fastas/*'

  output:
  path "*aligned.fasta"
  path "*.treefile"

  script:
  """
  # First concatenate all the fastas
  # How can I operate on all the fastas?
  cat fastas/*.fasta > concatenated.fasta
  mafft --auto concatenated.fasta > aligned.fasta
  iqtree -T AUTO -s aligned.fasta -B 1000
  """

}
