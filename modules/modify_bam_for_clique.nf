process MODIFY_BAM {

  publishDir "${params.outdir}/2_clique_snv/", mode:'copy', pattern:'*.{txt}'

  input:
  tuple val(sampleName), path ("${sampleName}.markdup.bam"), path ("${sampleName}.markdup.bam.bai")

  output:
  path ("*_modified.sam"), emit: MODIFY_BAM_out
  path "*.txt"

  script:
  """
  samtools stats ${sampleName}.markdup.bam > ${sampleName}.markdup.bam.stats
  SEQ=\$(grep 'reads mapped:' *.markdup.bam.stats | cut -f3)

  if [[ \$SEQ -gt 499 ]]
  then
    echo 'More than 500 reads (duplicates removed) mapped to best reference'> ${sampleName}_MODIFY_BAM_info.txt
    # convert to sam
    samtools view -h -O SAM -o ${sampleName}.sam ${sampleName}.markdup.bam

    # Remove @PG header line
    grep -v "@PG" ${sampleName}.sam > ${sampleName}_modified.sam

  else
    echo 'Less than 500 reads (duplicates removed) mapped to best reference'> ${sampleName}_MODIFY_BAM_info.txt
  fi
  """

}
