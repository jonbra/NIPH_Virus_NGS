nextflow.enable.dsl=2

nf_mod_path = "$baseDir/modules"
// Kraken viral db built on june 7 2022 downloaded from https://benlangmead.github.io/aws-indexes/k2.
// Custom set of HCV subtypes added. Same sequences as in HCVgenosubtypes_8.5.19_clean.fa
kraken_db = "$baseDir/Data/Kraken_db/"
kraken_hcv_db = "$baseDir/Data/Kraken_db/HCV/"
blast_db = "$baseDir/Data/Blast_db/HCVgenosubtypes_8.5.19_clean.fa"
ref_file = "$baseDir/Data/References/3a_D17763.fa"
params.cpus=4
params.outdir = params.outpath + "/results/"
params.taxid = 11102 // Hepacivirus default

include { FASTQC } from "$nf_mod_path/fastqc.nf"
include { TRIM } from "$nf_mod_path/trim.nf"
include { FASTQC as FASTQC_TRIM } from "$nf_mod_path/fastqc.nf"
include { KRAKEN2 } from "$nf_mod_path/kraken2.nf"
include { SUBSET_KRAKEN2 } from "$nf_mod_path/subset_kraken2.nf"
include { SPADES } from "$nf_mod_path/spades.nf"
include { MULTIQC } from "$nf_mod_path/multiqc.nf"
include { BLASTN } from "$nf_mod_path/blastn.nf"
//include { ABACAS } from "$nf_mod_path/abacas.nf"

workflow {
  reads = Channel
          .fromPath(params.samplelist)
          .splitCsv(header:true, sep:",")
          .map{ row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2))}
  
  FASTQC(reads, 'raw')
  TRIM(reads)
  FASTQC_TRIM(TRIM.out.TRIM_out, 'trimmed')
  KRAKEN2(TRIM.out.TRIM_out, kraken_db)

  // Run Spades if --skip_assembly is not active
  if (!params.skip_assembly) {
    SUBSET_KRAKEN2(TRIM.out.TRIM_out, KRAKEN2.out.report, KRAKEN2.out.classified_reads_assignment)
    SPADES(SUBSET_KRAKEN2.out.subset_reads_fastq)
    BLASTN(SPADES.out.scaffolds, blast_db)
    //ABACAS(SPADES.out.scaffolds, ref_file)
  }

  // MultiQC
  // Takes FastQC and Kraken output
  FILES_FOR_MULTIQC = FASTQC.out.FASTQC_out.collect { it[1] }.mix(
  FASTQC_TRIM.out.FASTQC_out.collect { it[1] }.mix(
    KRAKEN2.out.report)
  ).collect()
  MULTIQC(FILES_FOR_MULTIQC)
}
