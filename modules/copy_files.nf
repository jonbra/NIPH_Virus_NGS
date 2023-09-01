process COPY_FILES {

    container 'jonbra/viral_haplo:1.3'

    label 'tiny'

    publishDir "${params.outdir}/", mode: 'copy', pattern:'*.{csv}'

    input:
    path(samplelist)
    path(outdir)
    //path(parameter_file)

    script:
    """
    cp ${samplelist} "$baseDir/${outdir}"
    """
}