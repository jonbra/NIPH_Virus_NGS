nextflow.enable.dsl=2

pipeline_version = "dev"
params.pipeline_version = pipeline_version // For use in modules
params.cpu=4
nf_mod_path = "$baseDir/modules"

// **********************************************************************************

ref_file = "$baseDir/Data/References/HCVgenosubtypes_8.5.19_clean.fa"
//primer_bed = "$baseDir/util/swift_primers.bed"
//primer_master_file = "$baseDir/util/sarscov2_v2_masterfile.txt"
//nscTrim_primer_file = "$baseDir/util/swift_amplicon_nscTrim_750b.txt"

//params.ref_id = "NC_045512.2"
params.align_tool = "bowtie2"
params.outdir = params.outpath + "/results/"


//vars_under_obs_file = "$baseDir/util/variants.csv"
//params.check_variants_py = "check_variants_" + pipeline_version + ".py"
//params.plotting_py = "plotting_" + pipeline_version + ".py"

// **********************************************************************************


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
//include { REPORT } from "$nf_mod_path/report.nf"
include { TRIM } from "$nf_mod_path/trim.nf"
include { FASTQC as FASTQC_TRIM } from "$nf_mod_path/fastqc.nf"
include { INDEX } from "$nf_mod_path/index.nf"
include { BOWTIE } from "$nf_mod_path/bowtie2.nf"
include { TANOTI } from "$nf_mod_path/tanoti.nf"
include { CONSENSUS} from "$nf_mod_path/consensus.nf"
include { MODIFY_FASTA } from "$nf_mod_path/modify_fasta.nf"
include { MULTIQC } from "$nf_mod_path/multiqc.nf"
/*include { NSCTRIM } from "$nf_mod_path/nsctrim.nf"
include { FASTP } from "$nf_mod_path/fastp.nf"
include { FASTQC as FASTQC_CLEAN } from "$nf_mod_path/fastqc.nf"



include { BOWTIE2_INDEX; BOWTIE2_ALIGN } from "$nf_mod_path/bowtie2.nf"
include { BWA_INDEX; BWA_ALIGN } from "$nf_mod_path/bwa.nf"

include { PICARD_WGSMETRICS } from "$nf_mod_path/picard.nf"

include { SAMTOOLS_MPILEUP } from "$nf_mod_path/samtools.nf"

include { IVAR_VARIANTS; IVAR_CONSENSUS; CAT_CONSENSUS } from "$nf_mod_path/ivar.nf"
include { VARSCAN2_VARIANTS; VARSCAN2_CONSENSUS } from "$nf_mod_path/varscan2.nf"

include { PANGOLIN as PANGOLIN_IVAR } from "$nf_mod_path/lineage.nf"
include { NEXTCLADE as NEXTCLADE_IVAR } from "$nf_mod_path/lineage.nf"

include { NOISE_EXTRACTOR; FRAMESHIFT_FINDER } from "$nf_mod_path/tools.nf"
include { CHECK_VARIANTS } from "$nf_mod_path/checkvariants.nf"

include { GENERATE_REPORT; QC_PLOTS; NEXTCLADE_FOR_FHI } from "$nf_mod_path/reportgenerator.nf"
*/
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
            FASTQC_TRIM.out.FASTQC_out.collect { it[1] }
    ).collect()
    MULTIQC(FILES_FOR_MULTIQC)
}
    //NSCTRIM(reads, nscTrim_primer_file)
    //FASTP(NSCTRIM.out.NSCTRIM_out)
    //FASTQC_CLEAN(FASTP.out.FASTP_out, 'clean')
    /*
    if ( params.align_tool == "bowtie2") {
        BOWTIE2_INDEX(ref_file)
        BOWTIE2_ALIGN(FASTP.out.FASTP_out, ref_file, BOWTIE2_INDEX.out.BOWTIE2_INDEX_out)
        ALIGNED = BOWTIE2_ALIGN.out.BOWTIE2_ALIGN_out
    } else if ( params.align_tool == "bwa") {
        BWA_INDEX(ref_file)
        BWA_ALIGN(FASTP.out.FASTP_out, ref_file, BWA_INDEX.out.BWA_INDEX_out)
        ALIGNED = BWA_ALIGN.out.BWA_ALIGN_out
    }

    PICARD_WGSMETRICS(ALIGNED, ref_file)
    SAMTOOLS_MPILEUP(ALIGNED, ref_file)

    IVAR_VARIANTS(SAMTOOLS_MPILEUP.out.SAMTOOLS_MPILEUP_out, ref_file)
    IVAR_CONSENSUS(SAMTOOLS_MPILEUP.out.SAMTOOLS_MPILEUP_out, ref_file)
    PANGOLIN_IVAR(IVAR_CONSENSUS.out.IVAR_CONSENSUS_NREMOVED_out, 'ivar')
    NEXTCLADE_IVAR(IVAR_CONSENSUS.out.IVAR_CONSENSUS_NREMOVED_out, 'ivar')

    VARSCAN2_VARIANTS(SAMTOOLS_MPILEUP.out.SAMTOOLS_MPILEUP_out, ref_file)
    VARSCAN2_CONSENSUS(ALIGNED.join(VARSCAN2_VARIANTS.out.VARSCAN2_VARIANTS_out), ref_file)

    // Combine all consensus files into one file
    CAT_CONSENSUS(IVAR_CONSENSUS.out.IVAR_CONSENSUS_NREMOVED_out.collect { it[1] })

    // Other tools
    NOISE_EXTRACTOR(ALIGNED.collect { it[1..2] })
    FRAMESHIFT_FINDER(CAT_CONSENSUS.out.FASTA_out)

    // check_variant requires both the BAM and VCF files, so it will run at the end
    CHECK_VARIANTS(
        ALIGNED.collect { it[1..2] },
        IVAR_VARIANTS.out.IVAR_VARIANTS_out.collect { it[1..2] },
        vars_under_obs_file,
        Channel.fromPath(params.samplelist)
    )



    // Report generator and QC
    GENERATE_REPORT(
        Channel.fromPath(params.samplelist),
        Channel.fromPath("$params.outpath/pipeline_info.txt"),
        NSCTRIM.out.NSCTRIM_log.collect(),
        FASTP.out.FASTP_json.collect(),
        BOWTIE2_ALIGN.out.BOWTIE2_log.collect(),
        PICARD_WGSMETRICS.out.PICARD_WGSMETRICS_out.collect { it[1] },
        IVAR_CONSENSUS.out.IVAR_CONSENSUS_out.collect { it[1] },
        IVAR_VARIANTS.out.IVAR_BCFTOOLS_STATS_out.collect(),
        PANGOLIN_IVAR.out.PANGOLIN_out.collect { it[2] },
        NEXTCLADE_IVAR.out.NEXTCLADE_out.collect { it[2] }
        )

    QC_PLOTS(GENERATE_REPORT.out.GENERATE_REPORT_out)
    NEXTCLADE_FOR_FHI(NEXTCLADE_IVAR.out.NEXTCLADE_out.collect { it[2] })
    */
