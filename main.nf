nextflow.enable.dsl=2

kraken_db = "$baseDir/Data/Kraken_db/"
blast_db = "$baseDir/Data/Blast_db/HCVgenosubtypes_8.5.19_clean.fa"
ref_file = "$baseDir/Data/References/3a_D17763.fa"
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
include { SUBSET_KRAKEN2 } from "$nf_mod_path/subset_kraken2.nf"
include { SPADES } from "$nf_mod_path/spades.nf"
include { MULTIQC } from "$nf_mod_path/multiqc.nf"
include { BLASTN } from "$nf_mod_path/blastn.nf"
include { BLAST_PARSE } from "$nf_mod_path/blast_parse.nf"
//include { ABACAS } from "$nf_mod_path/abacas.nf"
include { INDEX } from "$nf_mod_path/index.nf"
include { DEDUP } from "$nf_mod_path/dedup.nf"
include { BOWTIE2 } from "$nf_mod_path/bowtie2.nf"
//include { TANOTI } from "$nf_mod_path/tanoti.nf"
//include { HCV_GLUE_SQL } from "$nf_mod_path/hcv_glue_sql.nf"

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
    SUBSET_KRAKEN2(TRIM.out.TRIM_out, KRAKEN2.out.report, KRAKEN2.out.classified_reads_assignment)
    SPADES(SUBSET_KRAKEN2.out.subset_reads_fastq)
    BLASTN(SPADES.out.scaffolds, blast_db)
    BLAST_PARSE(BLASTN.out.blastn_out, blast_db, SPADES.out.scaffolds)
    //ABACAS(SPADES.out.scaffolds, ref_file)
  }

  if (params.map_to_reference) {
    DEDUP(TRIM.out.TRIM_out)
    INDEX(blast_db)
    BOWTIE2(DEDUP.out.DEDUP_out, blast_db, INDEX.out.INDEX_out)
    //TANOTI(DEDUP.out.DEDUP_out, blast_db)
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
