process CLIQUE_SNV {

  publishDir "${params.outdir}/2_clique_snv/", mode:'copy', pattern:'*.{fasta,json,txt}'

  input:
  tuple val(sampleName), path ("${sampleName}.major.sorted.bam"), path ("${sampleName}.major.sorted.bam.bai")

  output:
  path ("*.fasta"), emit: CLIQUE_out
  path "*.json"
  path "*.txt"

  script:
  """
  samtools stats ${sampleName}.major.sorted.bam > ${sampleName}.major.sorted.bam.stats
  SEQ=\$(grep 'reads mapped:' ${sampleName}.major.sorted.bam.stats | cut -f3)

  if [[ \$SEQ -gt 999 ]]
  then
    # convert to sam
    samtools view -h -O SAM -o ${sampleName}.sam ${sampleName}.major.sorted.bam

    # Remove @PG header line
    grep -v "@PG" ${sampleName}.sam > ${sampleName}_modified.sam

    # Run CLIQUE_SNV
    # Det blir noe surr med output directory tror jeg. Det er noen mellomfiler den ikke finner...
    # Tror det handler om at den ikke kan lage en mappe i nextflow-kjøringen? Det fungerer å bruke outdir til /home/jonr...
    java -jar /usr/local/bin/clique-snv.jar -m snv-illumina -threads $task.cpus -outDir TEST -in ${sampleName}_modified.sam -fdf extended
  else
    echo 'Less than 1000 reads (duplicates removed) mapped to best reference'> ${sampleName}_CLIQUE_info.txt
  fi
  """

}
