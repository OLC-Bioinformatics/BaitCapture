process MERGE_MAPPING_RESULTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c85b516872f711516305474353432f10480f882d%3A93a0052cece81cca2539a50e43842b0e5aa614c7-0' :
        'biocontainers/mulled-v2-c85b516872f711516305474353432f10480f882d--93a0052cece81cca2539a50e43842b0e5aa614c7-0' }"

    input:
    tuple val(meta), path(aligncov), path(idxstats), path(kma_res)

    output:
    tuple val(meta), path('*.mapstats.tsv')        , emit: mapstats
    tuple val(meta), path('*.presence_absence.tsv'), emit: presence_absence
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kma = kma_res ? "--kma=${kma_res}" : ""
    """
    merge-mapping-results.R \
        --aligncov=${aligncov} \
        --idxstats=${idxstats} \
        --output_prefix=${prefix} \
        --outdir . \
        $kma \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n 1 | grep -Eo '[0-9.]+ ')
        r-readr: \$(grep 'readr' package-versions.txt | awk '{print \$3}')
        r-dplyr: \$(grep 'dplyr' package-versions.txt | awk '{print \$3}')
        r-tidyr: \$(grep 'tidyr' package-versions.txt | awk '{print \$3}')
        r-optparse: \$(grep 'optparse' package-versions.txt | awk '{print \$3}')
    END_VERSIONS
    """
}
