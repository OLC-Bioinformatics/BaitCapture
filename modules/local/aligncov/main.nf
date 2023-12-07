process ALIGNCOV {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aligncov:0.0.2--pyh7cba7a3_0' :
        'biocontainers/aligncov:0.0.2--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(sorted_bam)

    output:
    tuple val(meta), path("*_stats.tsv"), emit: stats
    tuple val(meta), path("*_depth.tsv") , emit: depth
    // path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    aligncov \\
        -i $sorted_bam \\
        -o ${prefix}
    """
}
