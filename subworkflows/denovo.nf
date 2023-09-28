//
// This sub-workflow maps the trimmed reads to the different HCV reference genomes,
// identifies the two references with most reads mapped to,
// and then re-maps the reads to the two references.
//


workflow DENOVO_WORKFLOW {

    take:
      blast_for_mapping
      blast_for_summarize
      agens

    main:

    MAP_TO_GENOTYPES(blast_for_mapping)

    // Plot the coverage of all the genotype mappings
    PLOT_COVERAGE(MAP_TO_GENOTYPES.out.DEPTH.collect())

    // Summarize the mapping statistics for all samples
    SUMMARIZE_MAPPING(MAP_TO_GENOTYPES.out.STATS.collect(),
                      MAP_TO_GENOTYPES.out.DEPTH.collect(),
                      blast_for_summarize,
                      agens)
}