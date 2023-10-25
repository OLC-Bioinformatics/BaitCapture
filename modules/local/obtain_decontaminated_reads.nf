process OBTAIN_DECONTAMINATED_READS {

    tag "OBTAIN_DECONTAMINATED_READS ${sample_id}"
    label 'process_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_1' :
        'quay.io/biocontainers/samtools:1.17--hd87286a_1' }"


    publishDir "${params.outdir}/decontaminated-reads/", pattern: "${sample_id}*.fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(host_aligned_sam)

    output:
    tuple val(sample_id), path("${sample_id}.host.aligned.sorted.bam"), emit: host_aligned_sorted_bam
    tuple val(sample_id), path("${sample_id}_R{1,2}.decontam.fastq"), emit: decontaminated_reads

    script:
    """
    # Convert SAM to sorted BAM
    samtools view --threads ${task.cpus} -b -f 0x0004 -f 0x0008 -f 0x0001 ${sample_id}.host.aligned.sam | samtools sort --threads ${task.cpus} > ${sample_id}.host.aligned.sorted.bam
    # Convert sorted BAM to FASTQ
    samtools fastq --threads ${task.cpus} -1 ${sample_id}_R1.decontam.fastq -2 ${sample_id}_R2.decontam.fastq -s shouldnotexist.fq ${sample_id}.host.aligned.sorted.bam
    """

}