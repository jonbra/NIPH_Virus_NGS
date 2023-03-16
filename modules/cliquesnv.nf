process CLIQUE_SNV {

  container 'jonbra/viral_haplo:1.3'

  publishDir "${params.outdir}/clique_snv/", mode:'copy', pattern:'*.{fasta,json}'

  label 'medium'

  input:
  path 'bamfiles/'

  output:
  path ("*.fasta"), emit: CLIQUE_out
  path "*.json"

  script:
  """
  #mkdir bamfiles/TMP
  cd bamfiles/
  for i in \$(ls *.bam)
  do
  BASE=\${i%_ref.fa.sorted.nodups.bam}.sam
  IN=\${i%_ref.fa.sorted.nodups.bam}_modified.sam
  
  # convert to sam
  samtools view -h -O SAM -o \${BASE} \${i}

  # Remove @PG header line
  grep -v "@PG" \${BASE} > \${IN}

  java -jar ${projectDir}/bin/clique-snv.jar -m snv-illumina -threads $task.cpus -outDir ../ -in \${IN} -fdf extended
  done
  """

}

//   #cliquesnv -m snv-illumina -threads $task.cpus -outDir ../TMP -in ${samfile} -fdf extended
