process MAP_MAJORITY_TANOTI {

    container 'jonbra/viral_haplo:1.3'

    label 'medium'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path references
    tuple val(sampleName), path ("${sampleName}.first_mapping.sorted.bam"), path ("${sampleName}.first_mapping.sorted.bam.bai")
    tuple val(sampleName), path(major_ref)

    publishDir "${params.outdir}/6_map", mode: 'copy', pattern:'*major*.{bam,bai}'
    publishDir "${params.outdir}/6_map", mode: 'copy', pattern:'*.{stats,log,sh,txt,yml,gz}'

    output:
    tuple val(sampleName), path ("${sampleName}.*.major.markdup.bam"), path ("${sampleName}.*.major.markdup.bam.bai"), optional: true, emit: majority_out
    tuple val(sampleName), path ("${sampleName}.*.major.markdup.bam")                                                , optional: true, emit: GLUE
    path "${sampleName}.*.major.markdup.bam_coverage.txt.gz"                                                         , optional: true, emit: DEPTH
    path "${sampleName}.*.major.markdup.bam.stats"                                                                   , optional: true, emit: STATS
    path "*.log"                                                                                                     , optional: true, emit: BOWTIE2_log
    path "*.{stats,sh,txt}"                                                                                          , optional: true

    script:
    """
    # Re-map against the reference with the most mapped reads
    major="\$(< ${major_ref})"

    #major="\$(samtools idxstats ${sampleName}.first_mapping.sorted.bam | cut -f 1,3 | sort -k2 -h | tail -1 | cut -f1)" # Reference with most reads
    #echo "\${major}" >> ${sampleName}.majority_mapping_MAPPING_info.txt
    samtools faidx ${references} "\${major}" > "\${major}".fa # Get the fasta sequence

    tanoti \
        -r "\${major}.fa" \
        -i ${read1} ${read2} \
        -o ${sampleName}.majority_mapping.sam \
        -p 1 -m ${params.tanoti_stringency_2}
    
    samtools sort -@ $task.cpus -T ${sampleName} -O bam -o ${sampleName}.majority_mapping.sorted.bam ${sampleName}.majority_mapping.sam
    samtools index ${sampleName}.majority_mapping.sorted.bam
    samtools stats ${sampleName}.majority_mapping.sorted.bam > ${sampleName}.majority_mapping.sorted.bam.stats

    SEQ2=\$(grep 'reads mapped:' ${sampleName}.majority_mapping.sorted.bam.stats| cut -f3)
    echo \${SEQ2} >> ${sampleName}.majority_mapping_MAPPING_info.txt

    # Then remove duplicates
    samtools sort -n ${sampleName}.majority_mapping.sorted.bam \
      | samtools fixmate -m - - \
      | samtools sort -O BAM \
      | samtools markdup --no-PG -r - ${sampleName}.\${major}.major.markdup.bam

    samtools index ${sampleName}.\${major}.major.markdup.bam

    samtools stats ${sampleName}.\${major}.major.markdup.bam > ${sampleName}.\${major}.major.markdup.bam.stats
    
    # Creating file with coverage per site for plotting later
    samtools depth -aa -d 1000000 ${sampleName}.\${major}.major.markdup.bam | gzip > ${sampleName}.\${major}.major.markdup.bam_coverage.txt.gz

    SEQ3=\$(grep 'reads mapped:' ${sampleName}.\${major}.major.markdup.bam.stats | cut -f3)
    echo \${SEQ3} >> ${sampleName}.majority_mapping_MAPPING_info.txt

    cp .command.sh ${sampleName}.tanoti.majority_mapping.sh

    cat <<-END_VERSIONS > tanoti.majority_mapping.versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}