process GLUE_PARSER {

    container 'jonbra/tidyverse_seqinr:2.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/8_summaries/", mode:'copy', pattern:'*.{tsv}'
    publishDir "${params.outdir}/logs/"       , mode:'copy', pattern:'*.{log,sh}'

    input:
    path (json)

    output:
    path '*tsv', optional:true, emit: GLUE_summary
    path '*{log,sh}'
    path "glue_parser_versions.yml", emit: versions

    script:
    """
    GLUE_json_parser.R ${json}

    cp .command.log glue_parser_command.log
    cp .command.sh glue_parser_command.sh

    cat <<-END_VERSIONS > glue_parser_versions.yml
    "${task.process}":
        R-session: \$(cat R_versions.txt)
    END_VERSIONS
    """


}