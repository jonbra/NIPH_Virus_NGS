process MAP {
    tag "$sampleName"
    errorStrategy 'ignore'

    label 'large'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path genome
    path genome_index

    publishDir "${params.outdir}/2_bam", mode: 'link', pattern:'*.{bam,bai}'
    publishDir "${params.outdir}/2_bam/log", mode: 'link', pattern:'*.{stats,log,sh}'

    output:
    tuple val(sampleName), path ("${sampleName}.sorted.bam"), path ("${sampleName}.sorted.bam.bai"), emit: MAP_out
    path "*.log", emit: BOWTIE2_log
    path "${sampleName}_QC_PASS"
    path "*.{stats,sh}"

    script:
    """
    bowtie2 \
        --threads $task.cpus \
        --local \
        --very-sensitive-local \
        -x $genome \
        -1 ${read1} -2 ${read2} \
        2> ${sampleName}.bowtie2.log \
        | samtools sort -@ $task.cpus -T ${sampleName} -O bam - > ${sampleName}.sorted.bam

    samtools index ${sampleName}.sorted.bam

    samtools stats ${sampleName}.sorted.bam > ${sampleName}.sorted.bam.stats

    SEQ=\$(grep 'raw total sequences:' ${sampleName}.sorted.bam.stats | cut -f3)
    if [ \$SEQ -gt 10 ] # 10000
    then
        echo 'PASS' > ${sampleName}_QC_PASS
    else
        echo 'FAIL' > ${sampleName}_QC_FAIL
    fi

    cp .command.sh ${sampleName}.bowtie2.sh

    # Get the reference with most mapped reads
    # Need to have the bai file in the same location
    # This gives the name of the most mapped reference.
    samtools idxstats Virus220220.sorted.bam | cut -f 1,3 | sort -nr -k 2 | head -1 | cut -f1

    # Then redo the mapping with only that.
    # One way is to pull it out of the genome reference file.
    # Another is to add all the individual genotypes into a genomes.config file like in nf-core.
    # In addition to the combined references. 
    """
}
/*
process MAP {
   publishDir "${params.outdir}/mapping", mode:'copy'

   input:
   path index
   tuple val(pair_id), path(trimmed_R)

   output:
   path "${pair_id}.sam"

   script:
   """
   bowtie2 --threads 8 -x index -S ${pair_id}.sam -1 ${name}_R1.fastq -2 ${name}_R2.fastq
   """
}
*/
