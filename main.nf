nextflow.enable.dsl=2

include { FASTQC }                from "./modules/fastqc.nf"
include { TRIM }                  from "./modules/trim.nf"
include { FASTQC as FASTQC_TRIM } from "./modules/fastqc.nf"
//include { KRAKEN2 }               from "./modules/kraken2.nf"
include { KRAKEN2_FOCUSED }       from "./modules/kraken2_focused.nf"
include { MULTIQC }               from "./modules/multiqc.nf"
include { SPADES }               from "./modules/spades.nf"
include { BLASTN }                from "./modules/blastn.nf"
include { BLAST_PARSE }           from "./modules/blast_parse.nf"
include { MAP_TO_GENOTYPES }      from "./modules/map_to_genotypes.nf"
include { PLOT_COVERAGE as PLOT_COVERAGE_MAJOR}         from "./modules/plot_coverage.nf"
include { PLOT_COVERAGE as PLOT_COVERAGE_MINOR }         from "./modules/plot_coverage.nf"
include { SUMMARIZE_MAPPING }     from "./modules/summarize_mapping.nf"
include { DEDUP }                 from "./modules/dedup.nf"
include { BOWTIE2 }               from "./modules/bowtie2.nf"
include { TANOTI }                from "./modules/tanoti.nf"
include { HCV_GLUE as HCV_GLUE_MAJOR}          from "./modules/hcv_glue.nf"
include { HCV_GLUE as HCV_GLUE_MINOR} from "./modules/hcv_glue.nf"
include { GLUE_PARSER }           from "./modules/glue_parser.nf"
include { MAP_ALL_REFERENCES}             from "./modules/map_all_references.nf"
include { MAP_MAJORITY_TANOTI }        from "./modules/map_majority_tanoti.nf"
include { MAP_MAJORITY_BOWTIE2 }        from "./modules/map_majority_bowtie2.nf"
include { MAP_MINORITY_TANOTI }        from "./modules/map_minority_tanoti.nf"
include { MAP_MINORITY_BOWTIE2 }        from "./modules/map_minority_bowtie2.nf"
include { IDENTIFY_MAJOR_MINOR }    from "./modules/identify_major_minor.nf"
//include { CONSENSUS_MAJOR }     from "./modules/consensus_major.nf"
//include { CONSENSUS_MINOR }     from "./modules/consensus_minor.nf"
include { SUMMARIZE_MAPPING }     from "./modules/summarize_mapping.nf"

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
      new_path = "${params.samplelist}".replaceAll('samplesheets/', '')
      file(params.samplelist).copyTo("${params.outdir}/${new_path}")
      //file(params.samplelist).copyTo("${params.outdir}/${params.samplelist}")
      // TODO: remove the characters "samplesheets" from the file path
      file(params.parameter_file).copyTo("${params.outdir}/${params.parameter_file}")
  }

  // Quality trim and taxonomic classification of the reads
  FASTQC(reads, 'raw')
  TRIM(reads)
  FASTQC_TRIM(TRIM.out.TRIM_out, 'trimmed')

  //KRAKEN2(TRIM.out.TRIM_out, params.kraken_all)
  KRAKEN2_FOCUSED(TRIM.out.TRIM_out, params.kraken_focused)

  // Run de novo assembly of the classified reads and identify the neares references in the Blast database
  SPADES(KRAKEN2_FOCUSED.out.classified_reads_fastq)
  BLASTN(SPADES.out.scaffolds, params.blast_db)
  BLAST_PARSE(BLASTN.out.blastn_out)

  // Map reads to all references using Bowtie2
  MAP_ALL_REFERENCES(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db) 

  // Summarize the mapping and identify major and minor subtypes
  IDENTIFY_MAJOR_MINOR(MAP_ALL_REFERENCES.out.idxstats, MAP_ALL_REFERENCES.out.DEPTH)

  // Map reads to the most abundant subtype, and a possible minor subtype
  if (params.mapper == "tanoti") {
      MAP_MAJORITY_TANOTI(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db, MAP_ALL_REFERENCES.out.sorted_out, IDENTIFY_MAJOR_MINOR.out.major_ref)
      MAP_MINORITY_TANOTI(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db, MAP_ALL_REFERENCES.out.sorted_out, IDENTIFY_MAJOR_MINOR.out.minor_ref)
      ch_glue_major = MAP_MAJORITY_TANOTI.out.GLUE
      ch_glue_minor = MAP_MINORITY_TANOTI.out.GLUE
      ch_depth_major = MAP_MAJORITY_TANOTI.out.DEPTH
      ch_depth_minor = MAP_MINORITY_TANOTI.out.DEPTH
  }
  else if (params.mapper == "bowtie2") {
      MAP_MAJORITY_BOWTIE2(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db, MAP_ALL_REFERENCES.out.sorted_out, IDENTIFY_MAJOR_MINOR.out.major_ref)
      MAP_MINORITY_BOWTIE2(KRAKEN2_FOCUSED.out.classified_reads_fastq, params.blast_db, MAP_ALL_REFERENCES.out.sorted_out, IDENTIFY_MAJOR_MINOR.out.minor_ref)
      ch_glue_major = MAP_MAJORITY_BOWTIE2.out.GLUE
      ch_glue_minor = MAP_MINORITY_BOWTIE2.out.GLUE
      ch_depth_major = MAP_MAJORITY_BOWTIE2.out.DEPTH
      ch_depth_minor = MAP_MINORITY_BOWTIE2.out.DEPTH
  }

  // Run GLUE analysis
  HCV_GLUE_MAJOR(ch_glue_major)
  HCV_GLUE_MINOR(ch_glue_minor)
  
  // Plot the coverage of all the genotype mappings
  PLOT_COVERAGE_MAJOR(ch_depth_major)
  PLOT_COVERAGE_MINOR(ch_depth_minor)

  // Summarize the mapping statistics for all samples
  if (params.mapper == "tanoti") {
    ch_stats = MAP_MAJORITY_TANOTI.out.STATS.collect().mix(MAP_MINORITY_TANOTI.out.STATS.collect())
    ch_depth = MAP_MAJORITY_TANOTI.out.DEPTH.collect().mix(MAP_MINORITY_TANOTI.out.DEPTH.collect())
  } else if (params.mapper == "bowtie2") {
    ch_stats = MAP_MAJORITY_BOWTIE2.out.STATS.collect().mix(MAP_MINORITY_BOWTIE2.out.STATS.collect())
    ch_depth = MAP_MAJORITY_BOWTIE2.out.DEPTH.collect().mix(MAP_MINORITY_BOWTIE2.out.DEPTH.collect())
  }

  // Summarize the mapping statistics for all samples
  SUMMARIZE_MAPPING(ch_stats,
                    ch_depth,
                    BLASTN.out.for_summarize.collect())


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
