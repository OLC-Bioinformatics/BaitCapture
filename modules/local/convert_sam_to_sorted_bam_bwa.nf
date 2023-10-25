process CONVERT_SAM_TO_SORTED_BAM_BWA {

    tag "CONVERT_SAM_TO_SORTED_BAM_BWA ${sample_id}"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_1' :
        'quay.io/biocontainers/samtools:1.17--hd87286a_1' }"

    publishDir "${params.outdir}/bwa/", pattern: "${sample_id}.bwa.aligned.sorted.bam", mode: 'copy'

    input:
    tuple val(sample_id), path(targets_aligned_sam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa.aligned.sorted.bam"), emit: aligned_sorted_bam

    script:
    """
    samtools view --threads ${task.cpus} -S -b ${sample_id}.bwa.aligned.sam > ${sample_id}.bwa.aligned.bam
    samtools sort --threads ${task.cpus} ${sample_id}.bwa.aligned.bam > ${sample_id}.bwa.aligned.sorted.bam
    """

}