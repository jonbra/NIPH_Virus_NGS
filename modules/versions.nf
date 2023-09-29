process VERSIONS {
  
  container 'ubuntu:22.04'

  input:
    path (versions)
  
  output:
    path "all_versions.txt"

  publishDir "${params.outdir}/versions", mode: 'copy', pattern:'all_versions.txt'
  
  script:
  """
  cat ${versions} > all_versions.txt
  """
}