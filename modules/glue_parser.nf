process GLUE_PARSER {

    container 'jonbra/tidyverse_seqinr:2.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/8_summaries/", mode:'copy', pattern:'*.{tsv}'
    publishDir "${params.outdir}/logs/"       , mode:'copy', pattern:'*.{log,sh}'

    input:
    path '*.json'

    output:
    path '*tsv', emit: GLUE_summary
    path '*{log,sh}'

    script:
    """
    GLUE_json_parser.R

    cp .command.log glue_parser_command.log
    cp .command.sh glue_parser_command.sh
    """


}