process HBV_RT_BLAST_PARSE {

    container 'jonbra/tidyverse_seqinr:1.0'

    errorStrategy 'terminate'

    label 'small'

    publishDir "${params.outdir}/9_HBV_resistance", mode:'copy', pattern:'*.tsv'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'

    input:
    path 'json_files/'
    path reference

    output:
    path "*.tsv"
    path "*.{log,sh}"

    script:
    """
    hbv_resistance.R ${reference}

    cp .command.log hbv_rt_blast_parse_command.log
    cp .command.sh hbv_rt_blast_parse_command.sh
    """
    
}