process MULTIQC {

    tag "MULTIQC ${multiqc_files}"
    label 'process_single'

    conda "bioconda::multiqc=1.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.15--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0' }"


    publishDir "${params.outdir}/multiqc/", mode: 'copy'

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data", emit: data
    path "*_plots", optional:true, emit: plots

    script:
    """
    multiqc -f ${projectDir}/assets/multiqc_config.yml ${multiqc_files}
    """

}