nextflow.enable.dsl=2

pipeline_version = "dev"
params.pipeline_version = pipeline_version // For use in modules
params.cpus=4
nf_mod_path = "$baseDir/modules"

// **********************************************************************************

ref_file = "$baseDir/Data/References/HCVgenosubtypes_8.5.19_clean.fa"

params.align_tool = "bowtie2"
params.outdir = params.outpath + "/results/"


File pipeline_tool_file = new File("$params.outpath/pipeline_info.txt")
pipeline_tool_file.write '\n' +
                         'Pipeline:\t' + pipeline_version + '\n' +
                         '\n' +
                         'RunFolder:\t' + params.outpath + '\n' +
                         'SampleSheet:\t' + params.samplelist + '\n' +
                         'Mapper:\t\t' + params.align_tool + '\n' +
                         '\n'


log.info """\
         HCV - N F   P I P E L I N E
         ===================================
         ref_file: $ref_file
         outpath : ${params.outpath}
         """
         .stripIndent()

include { FASTQC } from "$nf_mod_path/fastqc.nf"
include { TRIM } from "$nf_mod_path/trim.nf"
include { DEDUP } from "$nf_mod_path/dedup.nf"
include { FASTQC as FASTQC_TRIM } from "$nf_mod_path/fastqc.nf"
include { INDEX } from "$nf_mod_path/index.nf"
include { BOWTIE } from "$nf_mod_path/bowtie2.nf"
include { TANOTI } from "$nf_mod_path/tanoti.nf"
include { MODIFY_BAM } from "$nf_mod_path/modify_bam_for_clique.nf"
//include { CLIQUE_SNV } from "$nf_mod_path/cliquesnv.nf"
//include { CLUSTER } from "$nf_mod_path/clustering.nf"
include { CONSENSUS} from "$nf_mod_path/consensus.nf"
include { MODIFY_FASTA } from "$nf_mod_path/modify_fasta.nf"
include { MULTIQC } from "$nf_mod_path/multiqc.nf"
// Next add a blast step of the major consensus among all references. But need to collect all the consensus-sequences first

workflow {
    main:
    if (params.test) {
        reads = Channel
                .fromSRA(["SRR1762355", "SRR1762354"])
                .map{ tuple(it[0], it[1][0], it[1][1]) }
    }
    else {
        reads = Channel
                .fromPath(params.samplelist)
                .splitCsv(header:true, sep:",")
                .map{ row -> tuple(row.sample, file(row.fastq_1), file(row.fastq_2)) }
    }

    FASTQC(reads, 'raw')
    //REPORT(reads, 'raw')
    TRIM(reads)
    DEDUP(TRIM.out.TRIM_out)
    FASTQC_TRIM(DEDUP.out.DEDUP_out, 'trimmed')

    if ( params.align_tool == "bowtie2" ) {
        INDEX(ref_file)
        BOWTIE(DEDUP.out.DEDUP_out, ref_file, INDEX.out.INDEX_out)
        ALIGNED = BOWTIE.out.BOWTIE2_out
    } else if ( params.align_tool == "tanoti" ) {
        TANOTI(DEDUP.out.DEDUP_out, ref_file)
        ALIGNED = TANOTI.out.TANOTI_out
    }

    MODIFY_BAM(ALIGNED)
    //CLIQUE_SNV(MODIFY_BAM.out.MODIFY_BAM_out)
    //FILES_FOR_CLUSTER = CLIQUE_SNV.out.CLIQUE_out.collect()
    //CLUSTER(FILES_FOR_CLUSTER)
    CONSENSUS(ALIGNED, ref_file)
    MODIFY_FASTA(reads, CONSENSUS.out.CONSENSUS_fa)

    // MultiQC -- Needs input from all FastQC and fastp reports
    if ( params.align_tool == "bowtie2" ) {
      FILES_FOR_MULTIQC = FASTQC.out.FASTQC_out.collect { it[1] }.mix(
              FASTQC_TRIM.out.FASTQC_out.collect { it[1] }.mix(
                BOWTIE.out.BOWTIE2_log.collect { it[1] }
                )
      ).collect()
    } else if ( params.align_tool == "tanoti" ) {
      FILES_FOR_MULTIQC = FASTQC.out.FASTQC_out.collect { it[1] }.mix(
              FASTQC_TRIM.out.FASTQC_out.collect { it[1] }
      ).collect()
    }

    MULTIQC(FILES_FOR_MULTIQC)
}
