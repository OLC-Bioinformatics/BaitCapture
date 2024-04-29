process KMA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::kma=1.4.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.9--he4a0461_2' :
        'biocontainers/kma:1.4.9--he4a0461_2' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.sam"), emit: sam
    tuple val(meta), path("*.res"), emit: res
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    INDEX=`basename \$(find -L ./ -name "*.length.b") | sed 's/\\.length.b\$//'`
    [ -z "\$INDEX" ] && echo "KMA index files not found" 1>&2 && exit 1

    kma \\
        $args \\
        -sam \\
        -t ${task.cpus} \\
        -ipe ${reads} \\
        -t_db kma/\$INDEX \\
        -o ${prefix} \\
        > ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma -v))
    END_VERSIONS
    """
}
