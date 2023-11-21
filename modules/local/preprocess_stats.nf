process PREPROCESS_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::seqtk=1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_1' :
        'biocontainers/seqtk:1.4--he4a0461_1' }"

    input:
    tuple val(meta), path(raw_reads)
    tuple val(meta), path(preprocessed_reads)

    output:
    tuple val(meta), path("*_preprocess_stats.csv"), emit: preprocess_stats
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def forward_raw_reads = "${raw_reads[0]}"
    def reverse_raw_reads = "${raw_reads[1]}"
    def forward_preprocessed_reads = "${preprocessed_reads[0]}"
    def reverse_preprocessed_reads = "${preprocessed_reads[1]}"
    """
    pass_reads=`zcat ${forward_preprocessed_reads} | wc -l | awk '{print \$1}'`
    pass_reads=`expr \${pass_reads} / 2`
    total_reads=`zcat ${forward_raw_reads} | wc -l | awk '{print \$1}'`
    total_reads=`expr \${total_reads} / 2`
    prop_reads=`awk "BEGIN {print \${pass_reads}/\${total_reads}}"`

    pass_fwd_bases=`seqtk size ${forward_preprocessed_reads} | cut -f 2`
    pass_rev_bases=`seqtk size ${reverse_preprocessed_reads} | cut -f 2`
    pass_bases=`expr \${pass_fwd_bases} + \${pass_rev_bases}`
    total_fwd_bases=`seqtk size ${forward_raw_reads} | cut -f 2`
    total_rev_bases=`seqtk size ${reverse_raw_reads} | cut -f 2`
    total_bases=`expr \${total_fwd_bases} + \${total_rev_bases}`
    prop_bases=`awk "BEGIN {print \${pass_bases}/\${total_bases}}"`

    echo "sample,input_reads,reads_passing_preprocessing,proportion_reads_passing_preprocessing,input_bases,bases_passing_preprocessing,proportion_bases_passing_preprocessing" > ${prefix}_preprocess_stats.csv
    echo "${prefix},\${total_reads},\${pass_reads},\${prop_reads},\${total_bases},\${pass_bases},\${prop_bases}" >> ${prefix}_preprocess_stats.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}