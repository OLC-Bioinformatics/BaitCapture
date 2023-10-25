process BWA_INDEX_TARGETS {

    tag "BWA_INDEX_TARGETS ${targets}"
    label 'process_medium'

    conda "bioconda::bwa=0.7.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
        'quay.io/biocontainers/bwa:0.7.17--hed695b0_7' }"

    input:
    path targets

    output:
    tuple path(targets), path("*"), emit: indexed_targets

    script:
    """
    bwa index ${targets}
    """

}