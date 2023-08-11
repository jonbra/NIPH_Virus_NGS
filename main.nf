nextflow.enable.dsl=2

include { FASTQC }                from "./modules/fastqc.nf"
include { TRIM }                  from "./modules/trim.nf"
include { FASTQC as FASTQC_TRIM } from "./modules/fastqc.nf"
include { KRAKEN2 }               from "./modules/kraken2.nf"
include { KRAKEN2_FOCUSED }       from "./modules/kraken2_focused.nf"
include { SPADES }                from "./modules/spades.nf"
include { MULTIQC }               from "./modules/multiqc.nf"
include { BLASTN }                from "./modules/blastn.nf"
include { BLAST_PARSE }           from "./modules/blast_parse.nf"
include { MAP_TO_GENOTYPES }      from "./modules/map_to_genotypes.nf"
include { PLOT_COVERAGE }         from "./modules/plot_coverage.nf"
include { SUMMARIZE_MAPPING }     from "./modules/summarize_mapping.nf"
include { INDEX }                 from "./modules/index.nf"
include { DEDUP }                 from "./modules/dedup.nf"
include { BOWTIE2 }               from "./modules/bowtie2.nf"
include { HCV_GLUE_SQL }          from "./modules/hcv_glue.nf"
include { GLUE_PARSER }           from "./modules/glue_parser.nf"
include { CLIQUE_SNV }            from "./modules/cliquesnv.nf"

workflow {

  if (params.test) {
      reads = Channel
              .fromSRA('ERR10028751')
              .map{ tuple(it[0], it[1][0], it[1][1])}
  }
  else {
      reads = Channel
              .fromPath(params.samplelist)
              .splitCsv(header:true, sep:",")
              .map{ row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2))}
      
      // Make a copy of the samplelist and params-file in the result folder
      file(params.samplelist).copyTo("${params.outdir}/${params.samplelist}")
      //file(params.params-file).copyTo("${params.outdir}/${params.params-file}")
  }

  INDEX(params.blast_db) 
  FASTQC(reads, 'raw')
  
  TRIM(reads)
  FASTQC_TRIM(TRIM.out.TRIM_out, 'trimmed')
  KRAKEN2(TRIM.out.TRIM_out, params.kraken_all)
  KRAKEN2_FOCUSED(TRIM.out.TRIM_out, params.kraken_focused)

  // Only run SPADES and mapping if Kraken2 focused results in classified reads
  if(KRAKEN2_FOCUSED.out.classified_reads_fastq.exists()) {
    SPADES(KRAKEN2_FOCUSED.out.classified_reads_fastq)
    BLASTN(SPADES.out.scaffolds, params.blast_db)
    BLAST_PARSE(BLASTN.out.blastn_out)
    MAP_TO_GENOTYPES(BLAST_PARSE.out.FOR_MAPPING)
  }

  // Plot the coverage of all the genotype mappings
  PLOT_COVERAGE(MAP_TO_GENOTYPES.out.DEPTH.collect())

  // Summarize the mapping statistics for all samples
  SUMMARIZE_MAPPING(MAP_TO_GENOTYPES.out.STATS.collect(),
                    MAP_TO_GENOTYPES.out.DEPTH.collect(),
                    BLASTN.out.for_summarize.collect())

  // Run GLUE for HCV
  if (params.glue) {
    HCV_GLUE_SQL(MAP_TO_GENOTYPES.out.GLUE.collect())
    GLUE_PARSER(HCV_GLUE_SQL.out.GLUE_json)
  }

  // Run CliqueSNV
  if (params.cliquesnv) {
    CLIQUE_SNV(MAP_TO_GENOTYPES.out.GLUE.collect())
  }
  
  //
  // MultiQC
  //
  MULTIQC(
    params.multiqc_config,
    FASTQC.out.FASTQC_out.collect{it[1]}.ifEmpty([]),
    TRIM.out.log.collect{it[1]}.ifEmpty([]),
    KRAKEN2.out.report.collect{it[1]}.ifEmpty([]),
    KRAKEN2_FOCUSED.out.report.collect{it[1]}.ifEmpty([]),
    FASTQC_TRIM.out.FASTQC_out.collect{it[1]}.ifEmpty([]),
    MAP_TO_GENOTYPES.out.flagstat_dups.collect{it[1]}.ifEmpty([]),
    MAP_TO_GENOTYPES.out.flagstat_nodups.collect{it[1]}.ifEmpty([])
  )

}
