include { HBV_RT_BLAST }          from "../modules/hbv/hbv_rt_blast.nf"
include { HBV_RT_BLAST_PARSE }    from "../modules/hbv/hbv_rt_blast_parse.nf"

workflow HBV_WORKFLOW {

  take:


  main:
  
  
  emit:
    versions = ch_versions
}