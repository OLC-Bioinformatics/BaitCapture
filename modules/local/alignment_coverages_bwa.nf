process ALIGNMENT_COVERAGES_BWA {

    tag "ALIGNMENT_COVERAGES_BWA ${sample_id}"
    label 'process_single'

    conda "bioconda::aligncov=0.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aligncov:0.0.2--pyh7cba7a3_0' :
        'quay.io/biocontainers/aligncov:0.0.2--pyh7cba7a3_0' }"


    publishDir "${params.outdir}/bwa/", pattern: "${sample_id}*.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("*.tsv"), emit: alignment_coverages

    script:
    """
    aligncov -i ${sorted_bam} -o ${sample_id}.bwa
    """

}