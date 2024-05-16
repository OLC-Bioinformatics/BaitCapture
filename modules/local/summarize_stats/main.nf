process SUMMARIZE_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c85b516872f711516305474353432f10480f882d%3A93a0052cece81cca2539a50e43842b0e5aa614c7-0' :
        'biocontainers/mulled-v2-c85b516872f711516305474353432f10480f882d--93a0052cece81cca2539a50e43842b0e5aa614c7-0' }"

    input:
    tuple val(meta), path(raw_fastqscan_json), path(preprocessed_fastqscan_json), path(stats), path(fastp_log)

    output:
    tuple val(meta), path('*.sumstats.tsv')                 , emit: mapstats
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args               ?: ""
    def prefix       = task.ext.prefix             ?: "${meta.id}"
    def trimmed      = fastp_log                   ? "--trimmed=${fastp_log}" : ""
    def preprocessed = preprocessed_fastqscan_json ? "--preprocessed=${preprocessed_fastqscan_json}" : ""
    """
    summarize-statistics.R \
        --raw=${raw_fastqscan_json} \
        --stats=${stats} \
        --output_prefix=${prefix} \
        --outdir . \
        $trimmed \
        $preprocessed \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | grep -Eo '[0-9.]+ ')
        r-readr: \$(grep 'readr' package-versions.txt | awk '{print \$3}')
        r-dplyr: \$(grep 'dplyr' package-versions.txt | awk '{print \$3}')
        r-jsonlite: \$(grep 'jsonlite' package-versions.txt | awk '{print \$3}')
        r-stringr: \$(grep 'stringr' package-versions.txt | awk '{print \$3}')        
        r-tidyr: \$(grep 'tidyr' package-versions.txt | awk '{print \$3}')
        r-optparse: \$(grep 'optparse' package-versions.txt | awk '{print \$3}')
    END_VERSIONS
    """
}
