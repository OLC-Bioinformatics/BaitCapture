process TRIMMOMATIC {

    tag "TRIMMOMATIC ${sample_id}"
    label 'process_medium'

    conda "bioconda::trimmomatic=0.39"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.paired.trim*.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*.unpaired.trim_*.fastq.gz"), emit: unpaired_reads
    tuple val(sample_id), path("*.trimmomatic.log"), emit: log
    tuple val(sample_id), path("*.trimmomatic.summary"), emit: summary

    script:
    // def args = task.ext.args ?: ''
    // def qual_trim = task.ext.args2 ?: ''
    def prefix = "${sample_id}"
    def output = "${sample_id}.paired.trim_1.fastq.gz ${sample_id}.unpaired.trim_1.fastq.gz ${sample_id}.paired.trim_2.fastq.gz ${sample_id}.unpaired.trim_2.fastq.gz"
    def qual_trim = "${params.trimmomatic}"

    """
    trimmomatic PE \
        -threads ${task.cpus} \
        -summary ${sample_id}.trimmomatic.summary \
        ${reads} \
        ${output} \
        ${qual_trim} \
        2> ${sample_id}.trimmomatic.log
    """

}