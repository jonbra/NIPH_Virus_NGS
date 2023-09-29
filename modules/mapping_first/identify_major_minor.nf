process IDENTIFY_MAJOR_MINOR {
    
    container 'jonbra/tidyverse_seqinr:2.0'

    errorStrategy 'terminate'

    label 'small'

    input:
    tuple val(sampleName), path(idxstats)
    tuple val(sampleName), path(depth)

    output:
    tuple val(sampleName), path("${sampleName}.major_ref.txt"), optional: true, emit: major_ref
    tuple val(sampleName), path("${sampleName}.minor_ref.txt"), optional: true, emit: minor_ref
    path "identify_major_minor_versions.yml"                                  , emit: versions
    path "*{log,sh}"

    script:
    """
    summarize_mapping_to_all_references.R ${idxstats} ${depth} ${sampleName}

    cp .command.log identify_major_minor_command.log
    cp .command.sh identify_major_minor_command.sh

    cat <<-END_VERSIONS > identify_major_minor_versions.yml
    "${task.process}":
        R-session: \$(cat R_versions.txt)
    END_VERSIONS
    """

}