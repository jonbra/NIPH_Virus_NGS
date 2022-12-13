process MAP_TO_GENOTYPES {
    
    container 'jonbra/viral_haplo:1.2'

    label 'small'

    publishDir "${params.outdir}/6_map/", mode:'copy', pattern:'*.{bam,bai,yml,stats,log,gz}'

    input:
    // Her kan det være flere ref.fa-filer. Det kan variere
    tuple val(sampleName), path(references), path(read1), path(read2)

    output:
    tuple val(sampleName), path('*.bam'), path('*.bai'), emit: BAM
    path '*.stats'                                     , emit: STATS
    path '*.gz'                                        , emit: DEPTH
    path '*.yml'
    path '*.log'

    script:
    """
    # Index the genotype references
    for i in \$(ls *ref.fa) 
    do
    bowtie2-build \$i \$i
    
    bowtie2 \
        --threads $task.cpus \
        --local \
        --very-sensitive-local \
        -x \$i \
        -1 ${read1} -2 ${read2} \
        2> ${sampleName}.\$i.bowtie2.log \
        -S ${sampleName}.\$i.sam

    samtools sort -@ $task.cpus -O bam -o ${sampleName}.\$i.sorted.bam ${sampleName}.\$i.sam

    samtools index ${sampleName}.\$i.sorted.bam

    samtools stats ${sampleName}.\$i.sorted.bam | grep ^SN | cut -f 2- > ${sampleName}.\$i.sorted.bam.stats

    # Creating file with coverage per site for plotting later
    samtools depth -aa -d 1000000 ${sampleName}.\$i.sorted.bam | gzip > ${sampleName}.\$i.sorted.bam_coverage.txt.gz

    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}