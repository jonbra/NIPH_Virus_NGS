nextflow.enable.dsl=2

include { FASTQC }                from "./modules/fastqc.nf"
include { TRIM }                  from "./modules/trim.nf"
include { FASTQC as FASTQC_TRIM } from "./modules/fastqc.nf"
//include { KRAKEN2 }               from "./modules/kraken2.nf"
include { KRAKEN2_FOCUSED }       from "./modules/kraken2_focused.nf"
//include { INDEX }                 from "./modules/index.nf"
include { MULTIQC }               from "./modules/multiqc.nf"
include { MAP_TO_GENOTYPES }      from "./modules/map_to_genotypes.nf"
include { PLOT_COVERAGE }         from "./modules/plot_coverage.nf"
include { SUMMARIZE_MAPPING }     from "./modules/summarize_mapping.nf"
include { DEDUP }                 from "./modules/dedup.nf"
include { BOWTIE2 }               from "./modules/bowtie2.nf"
include { TANOTI }                from "./modules/tanoti.nf"
include { HCV_GLUE }          from "./modules/hcv_glue.nf"
include { GLUE_PARSER }           from "./modules/glue_parser.nf"
include { MAP_ALL}             from "./modules/map_all.nf"
//include { MAP_ALL_TANOTI }             from "./modules/map_all_tanoti.nf"
//include { MAP_ALL_BOWTIE2 }             from "./modules/map_all_bowtie2.nf"
include { MAP_MAJORITY_TANOTI }        from "./modules/map_majority_tanoti.nf"
include { MAP_MAJORITY_BOWTIE2 }        from "./modules/map_majority_bowtie2.nf"
include { MAP_MINORITY_TANOTI }        from "./modules/map_minority_tanoti.nf"
include { MAP_MINORITY_BOWTIE2 }        from "./modules/map_minority_bowtie2.nf"
//include { CONSENSUS_MAJOR }     from "./modules/consensus_major.nf"
//include { CONSENSUS_MINOR }     from "./modules/consensus_minor.nf"
//include { GLUE_REF_BASED }      from "./modules/glue_ref.nf"
//SUMMARIZE_REF_BASED

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

  // Quality trim and taxonomic classification of the reads
  FASTQC(reads, 'raw')
  TRIM(reads)
  FASTQC_TRIM(TRIM.out.TRIM_out, 'trimmed')

  //KRAKEN2(TRIM.out.TRIM_out, params.kraken_all)
  KRAKEN2_FOCUSED(TRIM.out.TRIM_out, params.kraken_focused)

  // Map reads
  MAP_ALL(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db) 
  if (params.mapper == "tanoti") {
      //MAP_ALL_TANOTI(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db) 
      MAP_MAJORITY_TANOTI(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db, MAP_ALL.out.sorted_out)
      MAP_MINORITY_TANOTI(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db, MAP_ALL.out.sorted_out)
      ch_glue = MAP_MAJORITY_TANOTI.out.GLUE.collect().mix(MAP_MINORITY_TANOTI.out.GLUE.collect())
      ch_depth = MAP_MAJORITY_TANOTI.out.DEPTH.collect().mix(MAP_MINORITY_TANOTI.out.DEPTH.collect())
      ch_stats = MAP_MAJORITY_TANOTI.out.STATS.collect().mix(MAP_MINORITY_TANOTI.out.STATS.collect())
      //  CONSENSUS_MAJOR(MAP_MAJORITY.out.majority_out)
      //  CONSENSUS_MINOR(MAP_MINORITY.out.minority_out)
  }
  else if (params.mapper == "bowtie2") {
      //MAP_ALL_BOWTIE2(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db, INDEX.out.INDEX_out) 
      MAP_MAJORITY_BOWTIE2(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db, INDEX.out.INDEX_out, MAP_ALL.out.sorted_out)
      MAP_MINORITY_BOWTIE2(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db, INDEX.out.INDEX_out, MAP_ALL.out.sorted_out)
      ch_glue = MAP_MAJORITY_BOWTIE2.out.GLUE.collect().mix(MAP_MINORITY_BOWTIE2.out.GLUE.collect())
      ch_depth = MAP_MAJORITY_BOWTIE2.out.DEPTH.collect().mix(MAP_MINORITY_BOWTIE2.out.DEPTH.collect())
      ch_stats = MAP_MAJORITY_BOWTIE2.out.STATS.collect().mix(MAP_MINORITY_BOWTIE2.out.STATS.collect())
      //  CONSENSUS_MAJOR(MAP_MAJORITY.out.majority_out)
      //  CONSENSUS_MINOR(MAP_MINORITY.out.minority_out)
  }


  HCV_GLUE(ch_glue)
  //SUMMARIZE_REF_BASED(.collect)
  
  // Plot the coverage of all the genotype mappings
  PLOT_COVERAGE(ch_depth)

  // Summarize the mapping statistics for all samples
  //SUMMARIZE_MAPPING(ch_stats,
  //                  ch_depth)
                    //BLASTN.out.for_summarize.collect())


  //
  // MultiQC
  //
  MULTIQC(
    params.multiqc_config,
    FASTQC.out.FASTQC_out.collect{it[1]}.ifEmpty([]),
    TRIM.out.log.collect{it[1]}.ifEmpty([]),
    //KRAKEN2.out.report.collect{it[1]}.ifEmpty([]),
    KRAKEN2_FOCUSED.out.report.collect{it[1]}.ifEmpty([]),
    FASTQC_TRIM.out.FASTQC_out.collect{it[1]}.ifEmpty([]),
    //MAP_TO_GENOTYPES.out.flagstat_dups.collect{it[1]}.ifEmpty([]),
    //MAP_TO_GENOTYPES.out.flagstat_nodups.collect{it[1]}.ifEmpty([])
  )

}
