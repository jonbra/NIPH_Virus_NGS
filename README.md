# learning_nextflow

>**Warning**
>This pipeline is still under active development and comes with no guarantees whatsoever. Please carefully inspect the output.

The purpose of this pipeline is to identify viral genotypes and, if possible, resistance mutations from Illumina data.  

## What does the pipeline do?
1. Taxonomic classification of reads.
2. De novo assembly of viral genomes.
3. Genotype identification based on blast search of assembled scaffolds against a user-defined database.
4. Map reads to the identified genotypes. 
5. Summarizes mapping statistics and genotype identification.
6. If available, checks for the presence of known resistance mutations.

## What does it not do?
1. Provide complete or accurate genome assemblies.
2. Provide haplotype-level assemblies.
3. Provide an accurate genotype quantification.

## How to run the pipeline
You need a Linux system with [Nextflow](https://www.nextflow.io/), [Docker](https://www.docker.com/) and [Git](https://git-scm.com/) installed and available in your `$PATH`. 

First clone the repo and enter the directory
```
git clone https://github.com/jonbra/learning_nextflow.git
cd learning_nextflow
```

Then you need to provide some parameters and store them in a file (e.g. `my_params.json`). Here's an example for HCV:
```
{
    "agens": "HCV",
    "//agens_comment": "",
    "samplelist": "samplesheets/2023.04.24_HCV_from_SRA.csv",
    "outdir": "HCV_SRA",
    "spades_mode": "rnaviral",
    "kraken_all": "\/media\/jonr\/SATA6TB\/Kraken_db\/", 
    "//kraken_all_comment": "PlusPFP database compiled by Ben Langmead on 9. December 2022",
    "kraken_focused": "\/home\/jonr\/Prosjekter\/learning_nextflow\/Data\/Kraken_db\/db_hepacivirus",  
    "//kraken_focused_comment": "NCBI Hepacivirus genomes with HCV genotypes from NCBI added",
    "blast_db": "\/home\/jonr\/Prosjekter\/learning_nextflow\/Data\/Blast_db\/HCVgenosubtypes_8.5.19_clean.fa", 
    "//blast_db_comment": "Custom set of HCV fasta files representing different genotypes. 2k1b recombinant added."
}
```

And then create a Samplesheet that lists the Illumina files to be analyzed. This has to be a comma-separated file (.csv) with three columns named: "sample", "fastq_1", and "fastq_2". The "sample" column contains the sample ids and the "fastq_1" and "fastq_2" columns contain the paths to the forward and reverse reads. Here's an example:
```
sample,fastq_1,fastq_2
Sample_01,/path/to/fastq/Sample_01_R1.fastq.gz,/path/to/fastq/Sample_01_R2.fastq.gz
Sample_02,/path/to/fastq/Sample_02_R1.fastq.gz,/path/to/fastq/Sample_02_R2.fastq.gz
```  

You can use the R script `create_samplesheet.R` in the bin directory. The fastq path can be to a directory with sub-directories of different fastq-samples:
```
Rscript bin/create_samplesheet.R /path/to/directory/with/sub-directories/with/fastq/files samplesheets/samplesheet.csv
```

Finally run the pipeline:
```
nextflow run main.nf --samplelist samplesheets/samplesheet.csv -params-file my_params.json
```

## Pipeline overview
