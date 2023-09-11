process MULTIQC {

    container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'

    label 'small'

    input:
    path 'multiqc_config.yaml'
    path ('fastqc_raw/*')
    path ('cutadapt/*')
    //path ('kraken2_all/*')
    path ('kraken2_subset/*')
    path ('fastqc_trimmed/*')
    path ('bowtie2_dups/*')
    path ('bowtie2_nodups/*')
    //path ('assembly_spades/*')

    output:
    path "*multiqc_report.html"     , emit: report
    path "*_data"                   , emit: data
    path 'multiqc_command.{log,sh}'
    path 'multiqc_versions.yml'

    publishDir "${params.outdir}/9_multiqc/", mode: 'copy', pattern:'*multiqc_*'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.yml'

    script:
    """
    multiqc . 
    mv multiqc_data 0_multiqc_data
    mv multiqc_report.html 0_multiqc_report.html

    cp .command.sh multiqc_command.sh
    cp .command.log multiqc_command.log

    cat <<-END_VERSIONS > multiqc_versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}


/* 
Viralrecon
process MULTIQC {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::multiqc=1.13a" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.13a--pyhdfd78af_1' :
        'quay.io/biocontainers/multiqc:1.13a--pyhdfd78af_1' }"

    input:
    path 'multiqc_config.yaml'
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path fail_reads_summary
    path fail_mapping_summary
    path 'amplicon_heatmap_mqc.tsv'
    path ('fastqc/*')
    path ('fastp/*')
    path ('kraken2/*')
    path ('bowtie2/*')
    path ('bowtie2/*')
    path ('ivar_trim/*')
    path ('picard_markduplicates/*')
    path ('mosdepth/*')
    path ('variants/*')
    path ('variants/*')
    path ('variants/*')
    path ('variants/*')
    path ('variants/*')
    path ('variants/*')
    path ('cutadapt/*')
    path ('assembly_spades/*')
    path ('assembly_unicycler/*')
    path ('assembly_minia/*')

    output:
    path "*multiqc_report.html"     , emit: report
    path "*_data"                   , emit: data
    path "*variants_metrics_mqc.csv", optional:true, emit: csv_variants
    path "*assembly_metrics_mqc.csv", optional:true, emit: csv_assembly
    path "*_plots"                  , optional:true, emit: plots
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom_config = multiqc_custom_config ? "--config $multiqc_custom_config" : ''
    """
    ## Run MultiQC once to parse tool logs
    multiqc -f $args $custom_config .

    ## Parse YAML files dumped by MultiQC to obtain metrics
    multiqc_to_custom_csv.py --platform illumina

    ## Manually remove files that we don't want in the report
    if grep -q ">skip_assembly<" workflow_summary_mqc.yaml; then
        rm -f *assembly_metrics_mqc.csv
    fi

    if grep -q ">skip_variants<" workflow_summary_mqc.yaml; then
        rm -f *variants_metrics_mqc.csv
    fi

    rm -f variants/report.tsv

    ## Run MultiQC a second time
    multiqc -f $args -e general_stats --ignore nextclade_clade_mqc.tsv $custom_config .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
*/