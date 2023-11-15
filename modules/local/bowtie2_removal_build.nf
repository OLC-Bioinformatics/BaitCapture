// Modifed from https://github.com/nf-core/mag/blob/ba7234959405b7a6df253d6d03158a7cde030c96/modules/local/bowtie2_removal_build.nf
process BOWTIE2_REMOVAL_BUILD {
    tag "$fasta"
    label 'process_medium'

    conda "bioconda::bowtie2=2.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.4.2--py38h1c8e9b9_1' :
        'biocontainers/bowtie2:2.4.2--py38h1c8e9b9_1' }"

    input:
    path fasta

    output:
    path 'bt2_index_base*', emit: index
    path "versions.yml"   , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bowtie
    bowtie2-build --threads $task.cpus $fasta "bt2_index_base"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
    END_VERSIONS
    """
}
