#!/usr/bin/env nextflow

/*
-------------------------------------------------------------------------------
    BaitCapture
-------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Define parameters
params.reads = ""
params.outdir = "${launchDir}/results"
params.ref = "${launchDir}/targets.fa"

// Create channels
reads_ch = Channel.fromFilePairs("${params.reads}", checkIfExists: true)
ref_ch = Channel.fromPath("${params.ref}", checkIfExists: true)

/*
-------------------------------------------------------------------------------
    PROCESSES
-------------------------------------------------------------------------------
*/

process FASTQC {

    cpus = 4

    tag "FASTQC ${reads}"

    publishDir "${params.outdir}/fastqc/", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip") , emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """

}

process TRIMMOMATIC {

    cpus = 4

    tag "TRIMMOMATIC ${reads}"

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
    def qual_trim = "ILLUMINACLIP:${projectDir}/assets/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

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

process MULTIQC {

    publishDir "${params.outdir}/multiqc/", mode: 'copy'

    input:
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    script:
    """
    multiqc -f ${projectDir}/assets/multiqc_config.yml ${multiqc_files}
    """

}



/*
-------------------------------------------------------------------------------
    WORKFLOW
-------------------------------------------------------------------------------
*/

workflow {

    // ref_ch.view()
    // reads_ch.view()

    FASTQC(reads_ch)

    TRIMMOMATIC(reads_ch)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files.collect().view()
    
    MULTIQC(ch_multiqc_files.collect())
   

}

// Print execution summary
workflow.onComplete {
   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}