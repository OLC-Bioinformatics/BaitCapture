process BWA_ALIGN_TARGETS {

    tag "BWA_ALIGN_TARGETS ${sample_id}"
    label 'process_medium'    

    conda "bioconda::bwa=0.7.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--he4a0461_11' :
        'quay.io/biocontainers/bwa:0.7.17--he4a0461_11' }"

    input:
    tuple path(targets), path("*"), val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bwa.aligned.sam"), emit: targets_aligned_sam

    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem -t ${task.cpus} \$INDEX ${reads} > ${sample_id}.bwa.aligned.sam
    """

}