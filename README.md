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

Then you need to provide some parameters...
```
```

...and create a Samplesheet that lists the Illumina files to be analyzed
```
```

## Pipeline overview
