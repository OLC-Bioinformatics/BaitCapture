process KMA_INDEX_TARGETS {

    tag "KMA_INDEX_TARGETS ${targets}"
    label 'process_medium'

    conda "bioconda::kma=1.4.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.9--he4a0461_2' :
        'quay.io/biocontainers/kma:1.4.9--he4a0461_2' }"

    input:
    path targets

    output:
    tuple path(targets), path("*"), emit: indexed_targets

    script:
    """
    kma index -i ${targets} -o targets.kma
    """

}