process DECONTAMINATION_STATS {

    tag "DECONTAMINATION_STATS ${sample_id}"
    label 'process_single'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_1' :
        'quay.io/biocontainers/samtools:1.17--hd87286a_1' }"


    publishDir "${params.outdir}/decontaminated-reads/", pattern: "${sample_id}_decontamination_stats.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_forward_reads), path(decontaminated_forward_reads)

    output:
    tuple val(sample_id), path("${sample_id}_decontamination_stats.csv"), emit: decontamination_stats

    script:
    """
    pass=`wc -l ${decontaminated_forward_reads} | awk '{print \$1}'`
    pass=`expr \${pass} / 2`
    total=`zcat ${trimmed_forward_reads} | wc -l | awk '{print \$1}'`
    total=`expr \${total} / 2`
    prop=`awk "BEGIN {print \${pass}/\${total}}"`
    echo "input_reads,reads_passing_contamination_filter,fraction_reads_passing_contamination_filter" > ${sample_id}_decontamination_stats.csv
    echo \${total},\${pass},\${prop} >> ${sample_id}_decontamination_stats.csv
    """

}