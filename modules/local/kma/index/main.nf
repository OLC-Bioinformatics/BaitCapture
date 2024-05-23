process KMA_INDEX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?                                                 
        'https://depot.galaxyproject.org/singularity/mulled-v2-4ddce7507400eb398d61467a4d4702dad2402183:f3f884ff6e32c46d776188595b3c395d867b600f-0' :
        'biocontainers/mulled-v2-4ddce7507400eb398d61467a4d4702dad2402183:f3f884ff6e32c46d776188595b3c395d867b600f-0' }"

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
        kma: \$(echo \$(kma -v) | sed 's/KMA-//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
