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
include { MAP_TO_GENOTYPES }      from "./modules/map_to_genotypes.nf"
include { PLOT_COVERAGE }         from "./modules/plot_coverage.nf"
include { SUMMARIZE_MAPPING }     from "./modules/summarize_mapping.nf"
include { INDEX }                 from "./modules/index.nf"
include { DEDUP }                 from "./modules/dedup.nf"
include { BOWTIE2 }               from "./modules/bowtie2.nf"
include { TANOTI }                from "./modules/tanoti.nf"
include { HCV_GLUE_SQL }          from "./modules/hcv_glue.nf"
include { GLUE_PARSER }           from "./modules/glue_parser.nf"
include { CLIQUE_SNV }            from "./modules/cliquesnv.nf"
include { HBV_RT_BLAST }          from "./modules/hbv_rt_blast.nf"
include { HBV_RT_BLAST_PARSE }    from "./modules/hbv_rt_blast_parse.nf"

// 
// Import sub-workflows
//
include { HCV_WORKFLOW }    from "./subworkflows/hcv.nf"
include { DENOVO_WORKFLOW } from "./subworkflows/denovo.nf"


// 
// Run main workflow
//
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
      file(params.parameter_file).copyTo("${params.outdir}/${params.parameter_file}")
  }

  INDEX(params.blast_db) 
  FASTQC(reads, 'raw')
  TRIM(reads)
  FASTQC_TRIM(TRIM.out.TRIM_out, 'trimmed')
  if (params.kraken_all) {
    KRAKEN2(TRIM.out.TRIM_out, params.kraken_all_db)
  }
  KRAKEN2_FOCUSED(TRIM.out.TRIM_out, params.kraken_focused)

  SPADES(KRAKEN2_FOCUSED.out.classified_reads_fastq)
  BLASTN(SPADES.out.scaffolds, params.blast_db)
  BLAST_PARSE(BLASTN.out.blastn_out)

  // Differentiate between HCV and other agents:
  if ( params.agens == 'HCV') {
    HCV_WORKFLOW(
      KRAKEN2_FOCUSED.out.classified_reads_fastq, 
      params.blast_db,
      BLASTN.out.for_summarize.collect()
      )
  } else {
    DENOVO_WORKFLOW(
      BLAST_PARSE.out.FOR_MAPPING,
      BLASTN.out.for_summarize.collect(),
      params.agens
    )
  }

  //
  // MultiQC
  //
  //MULTIQC(
  //  params.multiqc_config,
  //  FASTQC.out.FASTQC_out.collect{it[1]}.ifEmpty([]),
  //  TRIM.out.log.collect{it[1]}.ifEmpty([]),
  //  KRAKEN2.out.report.collect{it[1]}.ifEmpty([]),
  //  KRAKEN2_FOCUSED.out.report.collect{it[1]}.ifEmpty([]),
  //  FASTQC_TRIM.out.FASTQC_out.collect{it[1]}.ifEmpty([]),
  //  MAP_TO_GENOTYPES.out.flagstat_dups.collect{it[1]}.ifEmpty([]),
  //  MAP_TO_GENOTYPES.out.flagstat_nodups.collect{it[1]}.ifEmpty([])
  //)

}
/*


    //
    // MODULE: MultiQC
    //
    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowCommons.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config,
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_fail_reads_multiqc.ifEmpty([]),
            ch_fail_mapping_multiqc.ifEmpty([]),
            ch_amplicon_heatmap_multiqc.ifEmpty([]),
            FASTQC_FASTP.out.fastqc_raw_zip.collect{it[1]}.ifEmpty([]),
            FASTQC_FASTP.out.trim_json.collect{it[1]}.ifEmpty([]),
            ch_kraken2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_multiqc.collect{it[1]}.ifEmpty([]),
            ch_ivar_trim_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_flagstat_multiqc.collect{it[1]}.ifEmpty([]),
            ch_mosdepth_multiqc.collect{it[1]}.ifEmpty([]),
            ch_ivar_counts_multiqc.collect{it[1]}.ifEmpty([]),
            ch_bcftools_stats_multiqc.collect{it[1]}.ifEmpty([]),
            ch_snpeff_multiqc.collect{it[1]}.ifEmpty([]),
            ch_quast_multiqc.collect().ifEmpty([]),
            ch_pangolin_multiqc.collect{it[1]}.ifEmpty([]),
            ch_nextclade_multiqc.collect().ifEmpty([]),
            ch_cutadapt_multiqc.collect{it[1]}.ifEmpty([]),
            ch_spades_quast_multiqc.collect().ifEmpty([]),
            ch_unicycler_quast_multiqc.collect().ifEmpty([]),
            ch_minia_quast_multiqc.collect().ifEmpty([])
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
    */