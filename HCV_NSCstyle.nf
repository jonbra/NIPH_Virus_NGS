nextflow.enable.dsl=2

pipeline_version = "dev"
params.pipeline_version = pipeline_version // For use in modules
params.cpu=4
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
include { FASTQC as FASTQC_TRIM } from "$nf_mod_path/fastqc.nf"
include { INDEX } from "$nf_mod_path/index.nf"
include { BOWTIE } from "$nf_mod_path/bowtie2.nf"
include { TANOTI } from "$nf_mod_path/tanoti.nf"
include { CONSENSUS} from "$nf_mod_path/consensus.nf"
include { MODIFY_FASTA } from "$nf_mod_path/modify_fasta.nf"
include { MULTIQC } from "$nf_mod_path/multiqc.nf"

workflow {
    main:
    if (params.test) {
        reads = Channel
                .fromSRA(["SRR11939535", "SRR12473500"])
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
    FASTQC_TRIM(TRIM.out.CUT_out, 'trimmed')
    INDEX(ref_file)

    if ( params.align_tool == "bowtie2") {
        BOWTIE(TRIM.out.CUT_out, ref_file, INDEX.out.INDEX_out)
        ALIGNED = BOWTIE.out.BOWTIE2_out
    } else if ( params.align_tool == "tanoti") {
        TANOTI(TRIM.out.CUT_out, ref_file)
        ALIGNED = TANOTI.out.TANOTI_out
    }

    // Mapping mot entire database, mapping mot major og minor
    //CONSENSUS(MAP.out.MAP_out, ref_file)
    CONSENSUS(TRIM.out.CUT_out, ALIGNED, ref_file)
    MODIFY_FASTA(reads, CONSENSUS.out.CONSENSUS_fa)

    // MultiQC -- Needs input from all FastQC and fastp reports
    FILES_FOR_MULTIQC = FASTQC.out.FASTQC_out.collect { it[1] }.mix(
            FASTQC_TRIM.out.FASTQC_out.collect { it[1] }.mix(
              BOWTIE.out.BOWTIE2_log.collect { it[1] }
              )
    ).collect()
    MULTIQC(FILES_FOR_MULTIQC)
}
