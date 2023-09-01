process IDENTIFY_MAJOR_MINOR {
    
    container 'jonbra/tidyverse_seqinr:1.0'

    errorStrategy 'terminate'

    label 'small'

    input:
    tuple val(sampleName), path(idxstats)

    output:
    path "${sampleName}.major_ref.txt", optional: true, emit: major_ref
    path "${sampleName}.minor_ref.txt", optional: true, emit: minor_ref
    path "R_versions.txt"
    path "*{log,sh}"

    script:
    """
    summarize_mapping_to_all_references.R ${idxstats} ${sampleName}

    cp .command.log summarize_mapping_command.log
    cp .command.sh summarize_mapping_command.sh
    """

}