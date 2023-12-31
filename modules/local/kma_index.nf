process KMA_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "bioconda::kma=1.4.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.9--he4a0461_2' :
        'biocontainers/kma:1.4.9--he4a0461_2' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("kma"), emit: index
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir kma
    kma \\
        index \\
        $args \\
        -i $fasta \\
        -o kma/${fasta.baseName}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma -v))
    END_VERSIONS
    """
}
