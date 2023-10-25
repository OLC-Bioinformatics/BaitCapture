process KMA_ALIGN_TARGETS {

    tag "KMA_ALIGN_TARGETS ${sample_id}" 
    label 'process_medium'

    conda "bioconda::kma=1.4.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.9--he4a0461_2' :
        'quay.io/biocontainers/kma:1.4.9--he4a0461_2' }"

    input:
    tuple path(indexed_targets), path("*"), val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.kma.aligned.sam"), emit: targets_aligned_sam

    script:
    """
    kma -mem_mode -ex_mode -1t1 -vcf \
        -t ${task.cpus} \
        -ipe ${reads} \
        -t_db targets.kma \
        -o ${sample_id}.kma.aligned.temp \
        -sam \
        > ${sample_id}.kma.aligned.sam
    """

}