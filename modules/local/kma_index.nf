process KMA_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "bioconda::kma=1.4.14,bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4ddce7507400eb398d61467a4d4702dad2402183:338e3d4c5936de29a9734afd7b7ed69281705b1a-0' :
        'biocontainers/mulled-v2-4ddce7507400eb398d61467a4d4702dad2402183:338e3d4c5936de29a9734afd7b7ed69281705b1a-0' }"

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
