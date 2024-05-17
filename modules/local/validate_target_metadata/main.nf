process VALIDATE_TARGET_METADATA {
    label 'process_single'

    conda "bioconda::pandas=0.19.2 bioconda::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ec389033376a3033fa570b73c1569fb99993aaf4:7fb9f5a1996d57d6c45b8d0350ed82996a0e1071-0' :
        'biocontainers/mulled-v2-ec389033376a3033fa570b73c1569fb99993aaf4:7fb9f5a1996d57d6c45b8d0350ed82996a0e1071-0' }"  

    input:
    tuple val(meta), path(targets)
    path (target_metadata)

    output:
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    validate-target-metadata.py \\
        ${targets} \\
        ${target_metadata}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(grep "pandas" package-versions.txt | sed 's/pandas==//g')
        biopython: \$(grep "biopython" package-versions.txt | sed 's/biopython==//g')
    END_VERSIONS
    """
}