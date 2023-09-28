process MAP_MINORITY_TANOTI {

    container 'jonbra/viral_haplo:1.3'

    label 'medium'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path references
    tuple val(sampleName), path ("${sampleName}.first_mapping.sorted.bam"), path ("${sampleName}.first_mapping.sorted.bam.bai")
    tuple val(sampleName), path(minor_ref)

    publishDir "${params.outdir}/6_map", mode: 'copy', pattern:'*minor*.{bam,bai}'
    publishDir "${params.outdir}/6_map", mode: 'copy', pattern:'*.{stats,log,sh,txt,yml}'

    output:
    tuple val(sampleName), path ("${sampleName}.*.minor.markdup.bam"), path ("${sampleName}.*.minor.markdup.bam.bai"), optional: true, emit: minority_out
    tuple val(sampleName), path ("${sampleName}.*.minor.markdup.bam")                                                , optional: true, emit: GLUE
    path "${sampleName}.*.minor.markdup.bam_coverage.txt.gz"                                                         , optional: true, emit: DEPTH
    path "${sampleName}.*.minor.markdup.bam.stats"                                                                   , optional: true, emit: STATS
    path "*.log"                                                                                                     , optional: true, emit: BOWTIE2_log
    path "*.{stats,sh,txt}"                                                                                          , optional: true

    script:
    """
    # Map reads against the subtype with the second most mapped reads
    minor="\$(cat ${minor_ref} | head -1)" # Getting the subtype with the second most mapped reads
    reads="\$(cat ${minor_ref} | head -2 | tail -1)" # How many reads mapped to this subtype in the first mapping. After duplicate removal
    cov="\$(cat ${minor_ref} | head -3 | tail -1)" # The coverage of this subtype in the first mapping. After duplicate removal. Depth cutoff of 5 or more.

    # Write info to file
    echo "\${minor}" >> ${sampleName}.minority_mapping_MAPPING_info.txt
    echo "\${reads}" >> ${sampleName}.minority_mapping_MAPPING_info.txt
    echo "\${cov}" >> ${sampleName}.minority_mapping_MAPPING_info.txt

    # Map to the minor reference if it had more than minAgensRead reads and more than 5x coverage (min coverage threshold of 5) in the first mapping
    # Variables are wrapped in \$(()) to make them numbers and not characters
    if [ \$(("\${reads}")) -gt \$((${params.minAgensRead})) ] && [ \$(("\${cov}")) -gt \$((50)) ]; then 

    samtools faidx ${references} "\${minor}" > "\${minor}".fa # Get the fasta sequence of the minor reference

    tanoti \
        -r "\${minor}.fa" \
        -i ${read1} ${read2} \
        -o ${sampleName}.minority_mapping.sam \
        -p 1 -m ${params.tanoti_stringency_2}
    
    samtools sort -@ $task.cpus -T ${sampleName} -O bam -o ${sampleName}.minority_mapping.sorted.bam ${sampleName}.minority_mapping.sam
    samtools index ${sampleName}.minority_mapping.sorted.bam
    samtools stats ${sampleName}.minority_mapping.sorted.bam > ${sampleName}.minority_mapping.sorted.bam.stats

    SEQ2=\$(grep 'reads mapped:' ${sampleName}.minority_mapping.sorted.bam.stats | cut -f3)
    echo \${SEQ2} >> ${sampleName}.minority_mapping_MAPPING_info.txt

    # Then remove duplicates
    samtools sort -n ${sampleName}.minority_mapping.sorted.bam \
      | samtools fixmate -m - - \
      | samtools sort -O BAM \
      | samtools markdup --no-PG -r - ${sampleName}.\${minor}.minor.markdup.bam

    samtools index ${sampleName}.\${minor}.minor.markdup.bam

    samtools stats ${sampleName}.\${minor}.minor.markdup.bam > ${sampleName}.\${minor}.minor.markdup.bam.stats

    # Creating file with coverage per site for plotting later
    samtools depth -aa -d 1000000 ${sampleName}.\${minor}.minor.markdup.bam | gzip > ${sampleName}.\${minor}.minor.markdup.bam_coverage.txt.gz

    SEQ3=\$(grep 'reads mapped:' ${sampleName}.\${minor}.minor.markdup.bam | cut -f3)
    echo \${SEQ3} >> ${sampleName}.minority_mapping_MAPPING_info.txt

    fi

    cp .command.sh ${sampleName}.tanoti.minority_mapping.sh

    cat <<-END_VERSIONS > tanoti.minority_mapping.versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}