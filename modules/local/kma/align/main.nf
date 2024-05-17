process KMA_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?                                                 
        'https://depot.galaxyproject.org/singularity/mulled-v2-4ddce7507400eb398d61467a4d4702dad2402183%3Af3f884ff6e32c46d776188595b3c395d867b600f-0' :
        'biocontainers/mulled-v2-4ddce7507400eb398d61467a4d4702dad2402183%3Af3f884ff6e32c46d776188595b3c395d867b600f-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.res"), emit: res
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    """
    INDEX=`basename \$(find -L ./ -name "*.length.b") | sed 's/\\.length.b\$//'`
    [ -z "\$INDEX" ] && echo "KMA index files not found" 1>&2 && exit 1

    kma \\
        $args \\
        -sam \\
        -t ${task.cpus} \\
        -ipe ${reads} \\
        -t_db kma/\$INDEX \\
        -o ${prefix} \\
        | samtools $samtools_command $args2 -@ $task.cpus -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kma: \$(echo \$(kma -v) | sed 's/KMA-//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
