process MAP_TO_GENOTYPES {
    
    container 'jonbra/viral_haplo:1.3'

    label 'small'

    publishDir "${params.outdir}/6_map/", mode:'copy', pattern:'*.{bam,bai,stats,gz}'
    publishDir "${params.outdir}/logs/", mode:'copy', pattern:'*.{log,sh}'
    publishDir "${params.outdir}/versions/", mode:'copy', pattern:'*.yml'

    input:
    // Her kan det være flere ref.fa-filer. Det kan variere
    tuple val(sampleName), path(references), path(read1), path(read2)

    output:
    tuple val(sampleName), path("${sampleName}*nodups.bam"), path("${sampleName}*nodups.bam.bai"), emit: BAM
    tuple val(sampleName), path("*.withdups.flagstat")                                           , emit: flagstat_dups
    tuple val(sampleName), path("*.nodups.flagstat")                                             , emit: flagstat_nodups
    path '*.nodups.bam'                                                                          , emit: GLUE
    path '*.stats'                                                                               , emit: STATS
    path '*.gz'                                                                                  , emit: DEPTH
    path '*.yml'
    path '*.log'
    path '*.sh'

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

    # Create stats before removing duplicates
    samtools flagstat ${sampleName}.\$i.sorted.bam > ${sampleName}.\$i.sorted.withdups.flagstat

    # Remove duplicates
    samtools sort -n ${sampleName}.\$i.sorted.bam > ${sampleName}.\$i.sorted.byQuery.bam 
	samtools fixmate -m ${sampleName}.\$i.sorted.byQuery.bam ${sampleName}.\$i.sorted.byQuery.fix.bam
    samtools sort ${sampleName}.\$i.sorted.byQuery.fix.bam > ${sampleName}.\$i.sorted.byQuery.fix.sorted.bam
    samtools markdup -r ${sampleName}.\$i.sorted.byQuery.fix.sorted.bam ${sampleName}.\$i.sorted.nodups.bam
    samtools index ${sampleName}.\$i.sorted.nodups.bam

    # Create stats after removing duplicates
    samtools stats ${sampleName}.\$i.sorted.nodups.bam | grep ^SN | cut -f 2- > ${sampleName}.\$i.sorted.nodups.bam.stats
    samtools flagstat ${sampleName}.\$i.sorted.nodups.bam > ${sampleName}.\$i.sorted.nodups.flagstat

    # Creating file with coverage per site for plotting later
    samtools depth -aa -d 1000000 ${sampleName}.\$i.sorted.nodups.bam | gzip > ${sampleName}.\$i.sorted.nodups.bam_coverage.txt.gz

    done

    cp .command.log ${sampleName}.bowtie2_command.log
    cp .command.sh ${sampleName}.bowtie2_command.sh

    cat <<-END_VERSIONS > bowtie2_versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}