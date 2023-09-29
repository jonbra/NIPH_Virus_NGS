nextflow.enable.dsl=2

//
// Import modules to use in the main workflow
//
include { FASTQC }                from "./modules/fastqc.nf"
include { TRIM }                  from "./modules/trim.nf"
include { FASTQC as FASTQC_TRIM } from "./modules/fastqc.nf"
include { KRAKEN2 }               from "./modules/kraken2.nf"
include { KRAKEN2_FOCUSED }       from "./modules/kraken2_focused.nf"
include { SPADES }                from "./modules/spades.nf"
include { MULTIQC }               from "./modules/multiqc.nf"
include { BLASTN }                from "./modules/blastn.nf"
include { BLAST_PARSE }           from "./modules/blast_parse.nf"



include { INDEX }                 from "./modules/index.nf"
include { DEDUP }                 from "./modules/dedup.nf"
include { BOWTIE2 }               from "./modules/bowtie2.nf"
include { TANOTI }                from "./modules/tanoti.nf"
include { HCV_GLUE_SQL }          from "./modules/hcv_glue.nf"
include { GLUE_PARSER }           from "./modules/glue_parser.nf"
include { CLIQUE_SNV }            from "./modules/cliquesnv.nf"
include { HBV_RT_BLAST }          from "./modules/hbv_rt_blast.nf"
include { HBV_RT_BLAST_PARSE }    from "./modules/hbv_rt_blast_parse.nf"
include { VERSIONS }              from "./modules/versions.nf"

// 
// Import sub-workflows
//
include { HCV_WORKFLOW }    from "./subworkflows/hcv.nf"
include { DENOVO_WORKFLOW } from "./subworkflows/denovo.nf"


// 
// Run main workflow
//
workflow {

  ch_versions = Channel.empty()

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
      file(params.parameter_file).copyTo("${params.outdir}/${params.parameter_file}")
  }

  //
  // Run processes common to all agents
  //
  INDEX(params.blast_db)
  ch_versions = ch_versions.mix(INDEX.out.versions) 

  FASTQC(reads, 'raw')
  ch_versions = ch_versions.mix(FASTQC.out.versions) 

  TRIM(reads)
  ch_versions = ch_versions.mix(TRIM.out.versions) 

  FASTQC_TRIM(TRIM.out.TRIM_out, 'trimmed')
  if (params.kraken_all) {
    KRAKEN2(TRIM.out.TRIM_out, params.kraken_all_db)
  }
  KRAKEN2_FOCUSED(TRIM.out.TRIM_out, params.kraken_focused)
  ch_versions = ch_versions.mix(KRAKEN2_FOCUSED.out.versions) 

  SPADES(KRAKEN2_FOCUSED.out.classified_reads_fastq)
  ch_versions = ch_versions.mix(SPADES.out.versions) 

  BLASTN(SPADES.out.scaffolds, params.blast_db)
  ch_versions = ch_versions.mix(BLASTN.out.versions) 

  BLAST_PARSE(BLASTN.out.blastn_out)
  ch_versions = ch_versions.mix(BLAST_PARSE.out.versions) 

  //
  // Differentiate between HCV and other agents
  //
  if ( params.agens == 'HCV') {
    HCV_WORKFLOW(
      KRAKEN2_FOCUSED.out.classified_reads_fastq, 
      params.blast_db,
      BLASTN.out.for_summarize.collect()
    )
    ch_versions = ch_versions.mix(HCV_WORKFLOW.out.versions)
    // Create channel for MultiQC that can have the same name regardless of agens
    //ch_XX_multiqc = HCV_WORKFLOW.out.XX
  } else {
    DENOVO_WORKFLOW(
      BLAST_PARSE.out.FOR_MAPPING,
      BLASTN.out.for_summarize.collect(),
      params.agens
    )
    ch_versions = ch_versions.mix(DENOVO_WORKFLOW.out.versions)
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
    FASTQC_TRIM.out.FASTQC_out.collect{it[1]}.ifEmpty([])
    //MAP_TO_GENOTYPES.out.flagstat_dups.collect{it[1]}.ifEmpty([]),
    //MAP_TO_GENOTYPES.out.flagstat_nodups.collect{it[1]}.ifEmpty([])
  )
  ch_versions = ch_versions.mix(MULTIQC.out.versions)

  //
  // Report versions
  //
  ch_versions.view()
  VERSIONS(ch_versions.unique().collectFile(name: 'combined_versions.yml'))

}
