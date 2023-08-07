#!/usr/bin/env nextflow

/*
-----------------------------------------------------------------------------------------
    BaitCapture
-----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// Define parameters
params.reads = ""
params.outdir = "${launchDir}/results"
params.ref = "${launchDir}/targets.fa"

// Create channels
reads_ch = Channel.fromFilePairs("${params.reads}", checkIfExists: true)
ref_ch = Channel.fromPath("${params.ref}", checkIfExists: true)

// Processes
process FASTQC {

    tag "FASTQC ${reads}"

    publishDir "${params.outdir}/fastqc/html/", pattern: "*.html", mode: 'copy'
    publishDir "${params.outdir}/fastqc/zip/", pattern: "*.zip", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(reads), path("*.html"), path("*.zip")

    script:
    """
    fastqc ${reads}
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
    tuple val(sample_id), path("*.log"), emit: log
    tuple val(sample_id), path("*.summary"), emit: summary

    script:
    // def args = task.ext.args ?: ''
    def prefix = "${sample_id}"
    def output = "${sample_id}.paired.trim_1.fastq.gz ${sample_id}.unpaired.trim_1.fastq.gz ${sample_id}.paired.trim_2.fastq.gz ${sample_id}.unpaired.trim_2.fastq.gz"
    // def qual_trim = task.ext.args2 ?: ''
    def qual_trim = "ILLUMINACLIP:${projectDir}/assets/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

    """
    trimmomatic PE \
        -threads ${task.cpus} \
        -summary ${sample_id}.summary \
        -trimlog ${sample_id}.log \
        ${reads} \
        ${output} \
        ${qual_trim}
    """

}

// Workflow
workflow {

    // ref_ch.view()
    // reads_ch.view()

    FASTQC(reads_ch)
    TRIMMOMATIC(reads_ch)

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