nextflow.enable.dsl=2

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
include { HBV_RT_BLAST } from "./modules/hbv_rt_blast.nf"
include { HBV_RT_BLAST_PARSE } from "./modules/hbv_rt_blast_parse.nf"

workflow {
  reads = Channel
          .fromPath(params.samplelist)
          .splitCsv(header:true, sep:",")
          .map{ row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2))}
 
  //blast_file = Channel.fromPath( params.blast_db )

  FASTQC(reads, 'raw')
  TRIM(reads)
  FASTQC_TRIM(TRIM.out.TRIM_out, 'trimmed')
  KRAKEN2(TRIM.out.TRIM_out, params.kraken_all)
  KRAKEN2_FOCUSED(TRIM.out.TRIM_out, params.kraken_focused)

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

  // Run Genotyping for HBV
  if (params.agens == 'HBV') {
    HBV_RT_BLAST(scaffolds, params.rt_domain)
    HBV_RT_BLAST_PARSE()
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
  MAP_TO_GENOTYPES.out.flagstat.collect{it[1]}.ifEmpty([])
)

}/*


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