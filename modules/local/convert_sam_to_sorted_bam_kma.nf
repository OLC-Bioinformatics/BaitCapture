process CONVERT_SAM_TO_SORTED_BAM_KMA {

    tag "CONVERT_SAM_TO_SORTED_BAM_KMA ${sample_id}"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_1' :
        'quay.io/biocontainers/samtools:1.17--hd87286a_1' }"

    publishDir "${params.outdir}/kma/", pattern: "${sample_id}.kma.aligned.sorted.bam", mode: 'copy'

    input:
    tuple val(sample_id), path(targets_aligned_sam)

    output:
    tuple val(sample_id), path("${sample_id}.kma.aligned.sorted.bam"), emit: aligned_sorted_bam

    script:
    """
    samtools view --threads ${task.cpus} -S -b ${sample_id}.kma.aligned.sam > ${sample_id}.kma.aligned.bam
    samtools sort --threads ${task.cpus} ${sample_id}.kma.aligned.bam > ${sample_id}.kma.aligned.sorted.bam
    """

}