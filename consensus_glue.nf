nextflow.enable.dsl=2


include { CONSENSUS_MAJOR }     from "./subworkflows/reference-based/consensus_major.nf"
include { HCV_GLUE_SQL }          from "./modules/hcv_glue.nf"
include { GLUE_PARSER }          from "./modules/glue_parser.nf"


workflow {

  if (params.test) {
      reads = Channel
              .fromSRA('ERR10028751')
              .map{ tuple(it[0], it[1][0], it[1][1])}
  }
  else {
      bams = Channel
              .fromPath(params.samplelist)
              .splitCsv(header:true, sep:",")
              .map{ row -> file(row.value)}
  }
  
  //CONSENSUS_MAJOR(bams)
  HCV_GLUE_SQL(bams)
  //GLUE_PARSER(HCV_GLUE_SQL.out.GLUE_json.collect())


}
