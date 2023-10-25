process BWA_INDEX_HOST {

    tag "BWA_INDEX_HOST ${host}"
    label 'process_medium'

    conda "bioconda::bwa=0.7.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
        'quay.io/biocontainers/bwa:0.7.17--hed695b0_7' }"

    input:
    path host

    output:
    tuple path(host), path("*"), emit: indexed_host

    script:
    """
    bwa index ${host}
    """

}