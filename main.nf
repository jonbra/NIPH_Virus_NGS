nextflow.enable.dsl=2

kraken_db = params.kraken_db
blast_db = params.blast_db
ref_file = params.reference
genotypes = params.genotypes
nf_mod_path = "$baseDir/modules"
// Kraken viral db built on june 7 2022 downloaded from https://benlangmead.github.io/aws-indexes/k2.
// Custom set of HCV subtypes added. Same sequences as in HCVgenosubtypes_8.5.19_clean.fa

// Setting default values
params.skip_assembly = false
params.map_to_reference = false

include { FASTQC } from "$nf_mod_path/fastqc.nf"
include { TRIM } from "$nf_mod_path/trim.nf"
include { FASTQC as FASTQC_TRIM } from "$nf_mod_path/fastqc.nf"
include { KRAKEN2 } from "$nf_mod_path/kraken2.nf"
//include { SUBSET_KRAKEN2 } from "$nf_mod_path/subset_kraken2.nf"
include { REPAIR } from "$nf_mod_path/repair.nf"
include { SPADES } from "$nf_mod_path/spades.nf"
include { MULTIQC } from "$nf_mod_path/multiqc.nf"
include { BLASTN } from "$nf_mod_path/blastn.nf"
include { BLAST_PARSE } from "$nf_mod_path/blast_parse.nf"
//include { ABACAS } from "$nf_mod_path/abacas.nf"
include { INDEX } from "$nf_mod_path/index.nf"
include { DEDUP } from "$nf_mod_path/dedup.nf"
include { BOWTIE2 } from "$nf_mod_path/bowtie2.nf"
include { TANOTI } from "$nf_mod_path/tanoti.nf"
//include { HCV_GLUE_SQL } from "$nf_mod_path/hcv_glue_sql.nf"
include { RVA_GENO } from "$nf_mod_path/rva_genotyping.nf"

workflow {
  reads = Channel
          .fromPath(params.samplelist)
          .splitCsv(header:true, sep:",")
          .map{ row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2))}
  
  FASTQC(reads, 'raw')
  TRIM(reads)
  FASTQC_TRIM(TRIM.out.TRIM_out, 'trimmed')
  KRAKEN2(TRIM.out.TRIM_out, kraken_db)

  // Run Spades if --skip_assembly is not active
  if (!params.skip_assembly) {
    //SUBSET_KRAKEN2(TRIM.out.TRIM_out, KRAKEN2.out.report, KRAKEN2.out.classified_reads_assignment)
    //REPAIR(SUBSET_KRAKEN2.out.subset_reads_fastq)
    SPADES(KRAKEN2.out.classified_reads_fastq)
    BLASTN(SPADES.out.scaffolds, blast_db)
    
    // Combine all Spades scaffolds files
    //COMBINE_SCAFFOLDS(SPADES.out.scaffolds.collect { it[1] })
    // Combine all blast.out files
    //COMBINE_BLAST(BLASTN.out.blastn_out.collect { it[1] })

    Tanke: En mulighet er å vente med Blast_parse til alt er ferdig. Deretter collecte alle blast.out og alle scaffolds og deretter prosessere disse sammen i R_scriptet...
    BLAST_PARSE(BLASTN.out.blastn_out.collect { it[1] }, blast_db, SPADES.out.scaffolds.collect { it[1] })
    //ABACAS(BLAST_PARSE.out.subtypes_references, BLAST_PARSE.out.scaffolds_fasta, BLAST_PARSE.out.genotypes)
  }

  if (params.map_to_reference) {
    DEDUP(TRIM.out.TRIM_out)
    INDEX(ref_file)
    BOWTIE2(DEDUP.out.DEDUP_out, ref_file, INDEX.out.INDEX_out)
    TANOTI(DEDUP.out.DEDUP_out, ref_file)
  }

  if (params.genotype) {
    RVA_GENO(SPADES.out.scaffolds, genotypes)
  }

  //HCV_GLUE_SQL(SPADES.out.scaffolds)

  // MultiQC
  // Takes FastQC and Kraken output
  FILES_FOR_MULTIQC = FASTQC.out.FASTQC_out.collect { it[1] }.mix(
  FASTQC_TRIM.out.FASTQC_out.collect { it[1] }.mix(
    KRAKEN2.out.report)
  ).collect()
  MULTIQC(FILES_FOR_MULTIQC)
}
