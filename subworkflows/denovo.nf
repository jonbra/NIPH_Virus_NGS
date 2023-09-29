//
// This sub-workflow maps the trimmed reads to the various genotypes,
// identified in the blast of the spades scaffolds against the reference blast database.  
//

include { MAP_TO_GENOTYPES }  from "../modules/denovo/map_to_genotypes.nf"
include { PLOT_COVERAGE }     from "../modules/denovo/plot_coverage.nf"
include { SUMMARIZE_MAPPING } from "../modules/denovo/summarize_mapping.nf"

workflow DENOVO_WORKFLOW {

    take:
      blast_for_mapping
      blast_for_summarize
      agens

    main:

    ch_versions = Channel.empty()

    MAP_TO_GENOTYPES(blast_for_mapping)
    ch_versions = ch_versions.mix(MAP_TO_GENOTYPES.out.versions)

    // Plot the coverage of all the genotype mappings
    PLOT_COVERAGE(MAP_TO_GENOTYPES.out.DEPTH.collect())
    ch_versions = ch_versions.mix(PLOT_COVERAGE.out.versions)

    // Summarize the mapping statistics for all samples
    SUMMARIZE_MAPPING(MAP_TO_GENOTYPES.out.STATS.collect(),
                      MAP_TO_GENOTYPES.out.DEPTH.collect(),
                      blast_for_summarize,
                      agens)
    ch_versions = ch_versions.mix(SUMMARIZE_MAPPING.out.versions)

    emit:
      versions = ch_versions
}