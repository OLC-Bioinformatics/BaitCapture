process CAT_RESULTS {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path mapstats
    path presence_absence
    path presence_absence_clusters
    path sumstats

    output:
    path "mapstats.tsv"                 , emit: cat_mapstats
    path "presence_absence.tsv"         , emit: cat_presence_absence
    path "presence_absence_clusters.tsv", emit: cat_presence_absence_clusters, optional: true
    path "sumstats.tsv"                 , emit: cat_sumstats
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # If at least one file exists, create a new file with the header and concatenate all files
    # Array of file naming patterns
    patterns=("*.mapstats.tsv" "*.presence_absence.tsv" "*.presence_absence_clusters.tsv" "*.sumstats.tsv")

    # Iterate over the array
    for pattern in "\${patterns[@]}"; do
        # Concatenated file with the same naming pattern
        concat_file="\${pattern#*.}"

        # Check if at least one file with the current naming pattern exists
        if ls \$pattern 1> /dev/null 2>&1; then
            # If it does, get the header line from the first file
            head -n 1 \$(ls \$pattern | head -n 1) > "\$concat_file"

            # Skip the header line in all files and append them to the concatenated file
            for file in \$pattern; do
                tail -n +2 "\$file" >> "\$concat_file"
            done
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
