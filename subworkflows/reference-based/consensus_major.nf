process CONSENSUS_MAJOR {

    container 'andersenlabapps/ivar'

    input:
    tuple val(sampleName), path ("${sampleName}.major.markdup.bam"), path ("${sampleName}.major.markdup.bam.bai")

    publishDir "${params.outdir}/ref-based", mode: 'copy', pattern:'*.{fa}'
    publishDir "${params.outdir}/ref-based", mode: 'copy', pattern:'*.{sh,yml}'

    output:
    tuple val(sampleName), path ("${sampleName}*.fa"), optional: true
    path "*.{sh,yml}", optional: true

    script:
    """
    samtools mpileup -aa -A -d 10000 -Q 20 ${sampleName}.major.markdup.bam | ivar consensus -p ${sampleName}_major -m 5

    cp .command.sh ${sampleName}.consensus_major.sh

    cat <<-END_VERSIONS > ${sampleName}.consensus_major.versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        iVar: \$(echo \$(ivar version 2>&1))
    END_VERSIONS
    """
}
