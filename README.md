# learning_nextflow

Warning!: This is my attempt at learning Nextflow and at the same time trying to develop a pipeline for predicting HCV genotypes and haplotypes from Illumina data. This pipeline and the results should be treated with extreme caution.  

Issues:  
2022.04.08: The process CLIQUE_SNV does not produce an output and does not stop running. 

Example usage:
```

mkdir Bowtie_SRA
sudo nextflow run HCV_NSCstyle.nf --cpu 8 --use_docker --align_tool "bowtie2" --outpath "Bowtie_SRA" --test
```
