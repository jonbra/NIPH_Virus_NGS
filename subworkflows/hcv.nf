//
// This sub-workflow maps the trimmed reads to the different HCV reference genomes,
// identifies the two references with most reads mapped to,
// and then re-maps the reads to the two references.
//

include { MAP_ALL_REFERENCES }                   from "../modules/hcv/map_all_references.nf"
include { IDENTIFY_MAJOR_MINOR }                 from "../modules/hcv/identify_major_minor.nf"
include { MAP_MAJORITY_TANOTI }                  from "../modules/hcv/map_majority_tanoti.nf"
include { MAP_MAJORITY_BOWTIE2 }                 from "../modules/hcv/map_majority_bowtie2.nf"
include { MAP_MINORITY_TANOTI }                  from "../modules/hcv/map_minority_tanoti.nf"
include { MAP_MINORITY_BOWTIE2 }                 from "../modules/hcv/map_minority_bowtie2.nf"
include { HCV_GLUE as HCV_GLUE_MAJOR }           from "../modules/hcv/hcv_glue.nf"
include { HCV_GLUE as HCV_GLUE_MINOR }           from "../modules/hcv/hcv_glue.nf"
include { GLUE_PARSER as HCV_GLUE_PARSER_MAJOR } from "../modules/hcv/glue_parser.nf" 
include { GLUE_PARSER as HCV_GLUE_PARSER_MINOR } from "../modules/hcv/glue_parser.nf" 
include { PLOT_COVERAGE as PLOT_COVERAGE_MAJOR } from "../modules/hcv/plot_coverage.nf"
include { PLOT_COVERAGE as PLOT_COVERAGE_MINOR } from "../modules/hcv/plot_coverage.nf"
include { SUMMARIZE_MAPPING_HCV }                from "../modules/hcv/summarize_mapping_hcv.nf"
//include { CLIQUE_SNV }                           from "../modules/hcv/clique_snv.nf"

workflow HCV_WORKFLOW {

  take:
    classified_reads
    blast_db
    blast_out

  main:
  
  ch_versions = Channel.empty()

  // Map reads to all references using Bowtie2
  MAP_ALL_REFERENCES(classified_reads, blast_db) 
  ch_versions = ch_versions.mix(MAP_ALL_REFERENCES.out.versions)

  // Summarize the mapping and identify major and minor subtypes
  IDENTIFY_MAJOR_MINOR(MAP_ALL_REFERENCES.out.idxstats, MAP_ALL_REFERENCES.out.DEPTH)
  ch_versions = ch_versions.mix(IDENTIFY_MAJOR_MINOR.out.versions)

  // Map reads to the most abundant subtype, and a possible minor subtype
  if (params.mapper == "tanoti") {
      MAP_MAJORITY_TANOTI(classified_reads, blast_db, MAP_ALL_REFERENCES.out.sorted_out, IDENTIFY_MAJOR_MINOR.out.major_ref)
      ch_versions = ch_versions.mix(MAP_MAJORITY_TANOTI.out.versions)
      MAP_MINORITY_TANOTI(classified_reads, blast_db, MAP_ALL_REFERENCES.out.sorted_out, IDENTIFY_MAJOR_MINOR.out.minor_ref)
      ch_glue_major = MAP_MAJORITY_TANOTI.out.GLUE
      ch_glue_minor = MAP_MINORITY_TANOTI.out.GLUE
      ch_depth_major = MAP_MAJORITY_TANOTI.out.DEPTH
      ch_depth_minor = MAP_MINORITY_TANOTI.out.DEPTH
  }
  else if (params.mapper == "bowtie2") {
      MAP_MAJORITY_BOWTIE2(classified_reads, blast_db, MAP_ALL_REFERENCES.out.sorted_out, IDENTIFY_MAJOR_MINOR.out.major_ref)
      ch_versions = ch_versions.mix(MAP_MAJORITY_BOWTIE2.out.versions)
      MAP_MINORITY_BOWTIE2(classified_reads, blast_db, MAP_ALL_REFERENCES.out.sorted_out, IDENTIFY_MAJOR_MINOR.out.minor_ref)
      ch_glue_major = MAP_MAJORITY_BOWTIE2.out.GLUE
      ch_glue_minor = MAP_MINORITY_BOWTIE2.out.GLUE
      ch_depth_major = MAP_MAJORITY_BOWTIE2.out.DEPTH
      ch_depth_minor = MAP_MINORITY_BOWTIE2.out.DEPTH
  }

  // Run GLUE analysis
  HCV_GLUE_MAJOR(ch_glue_major)
  HCV_GLUE_PARSER_MAJOR(HCV_GLUE_MAJOR.out.GLUE_json)
  ch_versions = ch_versions.mix(HCV_GLUE_PARSER_MAJOR.out.versions)
  HCV_GLUE_MINOR(ch_glue_minor)
  HCV_GLUE_PARSER_MINOR(HCV_GLUE_MINOR.out.GLUE_json)
  
  // Plot the coverage of all the genotype mappings
  PLOT_COVERAGE_MAJOR(ch_depth_major)
  ch_versions = ch_versions.mix(PLOT_COVERAGE_MAJOR.out.versions)
  PLOT_COVERAGE_MINOR(ch_depth_minor)

  // Summarize the mapping statistics for all samples
  if (params.mapper == "tanoti") {
    ch_stats = MAP_MAJORITY_TANOTI.out.STATS.collect().mix(MAP_MINORITY_TANOTI.out.STATS.collect())
    ch_depth = MAP_MAJORITY_TANOTI.out.DEPTH.collect().mix(MAP_MINORITY_TANOTI.out.DEPTH.collect())
  } else if (params.mapper == "bowtie2") {
    ch_stats = MAP_MAJORITY_BOWTIE2.out.STATS.collect().mix(MAP_MINORITY_BOWTIE2.out.STATS.collect())
    ch_depth = MAP_MAJORITY_BOWTIE2.out.DEPTH.collect().mix(MAP_MINORITY_BOWTIE2.out.DEPTH.collect())
  }

  SUMMARIZE_MAPPING_HCV(ch_stats,
                        ch_depth,
                        blast_out,
                        HCV_GLUE_PARSER_MAJOR.out.GLUE_summary.collect().mix(HCV_GLUE_PARSER_MINOR.out.GLUE_summary.collect()))
  ch_versions = ch_versions.mix(SUMMARIZE_MAPPING_HCV.out.versions)

  // Run CliqueSNV
/*   if (params.cliquesnv) {
    CLIQUE_SNV(MAP_TO_GENOTYPES.out.GLUE.collect())
  } */

  emit:
    versions = ch_versions
}