nextflow.enable.dsl=2

kraken_main = params.kraken_all
kraken_sub = params.kraken_focused
blast_file = params.blast_db

include { FASTQC } from "./modules/fastqc.nf"
include { TRIM } from "./modules/trim.nf"
include { FASTQC as FASTQC_TRIM } from "./modules/fastqc.nf"
include { KRAKEN2 } from "./modules/kraken2.nf"
include { KRAKEN2_FOCUSED } from "./modules/kraken2_focused.nf"
include { REPAIR } from "./modules/repair.nf"
include { SPADES } from "./modules/spades.nf"
include { MULTIQC } from "./modules/multiqc.nf"
include { BLASTN } from "./modules/blastn.nf"
include { BLAST_PARSE } from "./modules/blast_parse.nf"
include { MAP_TO_GENOTYPES } from "./modules/map_to_genotypes.nf"
include { PLOT_COVERAGE } from "./modules/plot_coverage.nf"
include { SUMMARIZE_MAPPING } from "./modules/summarize_mapping.nf"
include { ABACAS } from "./modules/abacas.nf"
include { INDEX } from "./modules/index.nf"
include { DEDUP } from "./modules/dedup.nf"
include { BOWTIE2 } from "./modules/bowtie2.nf"
include { TANOTI } from "./modules/tanoti.nf"
//include { HCV_GLUE_SQL } from "./modules/hcv_glue_sql.nf"
include { RVA_GENO } from "./modules/rva_genotyping.nf"

workflow {
  reads = Channel
          .fromPath(params.samplelist)
          .splitCsv(header:true, sep:",")
          .map{ row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2))}
 
  //blast_file = Channel.fromPath( params.blast_db )

  FASTQC(reads, 'raw')
  TRIM(reads)
  FASTQC_TRIM(TRIM.out.TRIM_out, 'trimmed')
  //KRAKEN2(TRIM.out.TRIM_out, kraken_main)
  KRAKEN2_FOCUSED(TRIM.out.TRIM_out, kraken_sub)

  // Run Spades if --skip_assembly is not active
  if (!params.skip_assembly) {
    //SUBSET_KRAKEN2(TRIM.out.TRIM_out, KRAKEN2.out.report, KRAKEN2.out.classified_reads_assignment)
    //REPAIR(SUBSET_KRAKEN2.out.subset_reads_fastq)
    SPADES(KRAKEN2_FOCUSED.out.classified_reads_fastq)
    BLASTN(SPADES.out.scaffolds, params.blast_db)
    BLAST_PARSE(BLASTN.out.blastn_out)
    MAP_TO_GENOTYPES(BLAST_PARSE.out.FOR_MAPPING)
    ABACAS(BLAST_PARSE.out.FOR_ABACAS)

    // Plot the coverage of all the genotype mappings
    PLOT_COVERAGE(MAP_TO_GENOTYPES.out.DEPTH.collect())

    // Summarize the mapping statistics for all samples
    SUMMARIZE_MAPPING(MAP_TO_GENOTYPES.out.STATS.collect())
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
    KRAKEN2_FOCUSED.out.report)
  ).collect()
  MULTIQC(FILES_FOR_MULTIQC)
}
