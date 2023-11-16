process DECONTAMINATION_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--hd87286a_0' :
        'biocontainers/samtools:1.18--hd87286a_0' }"

    input:
    tuple val(meta), path(forward_reads), path(decontaminated_forward_reads)

    output:
    tuple val(meta), path("*_decontamination_stats.csv"), emit: decontamination_stats
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pass=`zcat ${decontaminated_forward_reads} | wc -l | awk '{print \$1}'`
    pass=`expr \${pass} / 2`
    total=`zcat ${forward_reads} | wc -l | awk '{print \$1}'`
    total=`expr \${total} / 2`
    prop=`awk "BEGIN {print \${pass}/\${total}}"`
    echo "input_reads,reads_passing_contamination_filter,fraction_reads_passing_contamination_filter" > ${prefix}_decontamination_stats.csv
    echo \${total},\${pass},\${prop} >> ${prefix}_decontamination_stats.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}