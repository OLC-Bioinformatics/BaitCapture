#!/usr/bin/env nextflow

/*
-------------------------------------------------------------------------------
    BaitCapture
-------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include {validateParameters; paramsHelp; paramsSummaryLog} from 'plugin/nf-validation'

// Define parameters
params.reads = ""
params.outdir = "${launchDir}/results"
params.targets = "${launchDir}/targets.fa"
params.trimmomatic = "ILLUMINACLIP:${projectDir}/assets/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run main.nf --reads \"*_R{1,2}_001.fastq.gz\" --targets targets.fa")
   exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// Create channels
reads_ch = Channel.fromFilePairs("${params.reads}", checkIfExists: true)
targets_ch = Channel.fromPath("${params.targets}", checkIfExists: true)

/*
-------------------------------------------------------------------------------
    PROCESSES
-------------------------------------------------------------------------------
*/

process FASTQC {

    conda "bioconda::fastqc=0.11.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'biocontainers/fastqc:0.11.9--0' }"

    cpus = 4
    tag "FASTQC ${sample_id}"
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

    conda "bioconda::trimmomatic=0.39"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2':
        'biocontainers/trimmomatic:0.39--hdfd78af_2' }"

    cpus = 4
    tag "TRIMMOMATIC ${sample_id}"

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

process MULTIQC {

    withName: MULTIQC {
        conda "bioconda::multiqc=1.15"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/multiqc:1.15--pyhdfd78af_0' :
            'biocontainers/multiqc:1.15--pyhdfd78af_0' }"

    tag "MULTIQC ${multiqc_files}"
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

process BWA_INDEX {

    withName: BWA_INDEX {
        conda "bioconda::bwa=0.7.17"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
            'biocontainers/bwa:0.7.17--hed695b0_7' }"

    tag "BWA_INDEX ${targets}"

    input:
    path targets

    output:
    tuple path(targets), path("*"), emit: indexed_targets

    script:
    """
    bwa index ${targets}
    """

}

process BWA_ALIGN {

    withName BWA_ALIGN {
        conda "bioconda::bwa=0.7.17"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/bwa:0.7.17--h5bf99c6_8' :
            'biocontainers/bwa:0.7.17--h5bf99c6_8' }"

    cpus = 8
    tag "BWA_ALIGN ${sample_id}"
    publishDir "${params.outdir}/bwa-alignments/", pattern: "${sample_id}.aligned.sorted.bam", mode: 'copy'

    input:
    tuple path(targets), path("*"), val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.aligned.sorted.bam"), emit: aligned_targets

    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem -t ${task.cpus} \$INDEX ${reads} > ${sample_id}.aligned.sam
    samtools view --threads ${task.cpus} -S -b ${sample_id}.aligned.sam > ${sample_id}.aligned.bam
    samtools sort --threads ${task.cpus} ${sample_id}.aligned.bam > ${sample_id}.aligned.sorted.bam
    """

}

// process ALIGNMENT_STATS {

//     tag "ALIGNMENT_STATS ${sample_id}"
//     publishDir "${params.outdir}/alignment-stats/", pattern: "${sample_id}.sorted.bam.stats.tsv", mode: 'copy'

//     input:
//     tuple val(sample_id), path(sorted_bam)

//     output:
//     tuple val(sample_id), path("${sample_id}.sorted.bam.stats.tsv"), emit: alignment_stats

//     script:
//     """
//     samtools index ${sorted_bam}
//     samtools stats ${sorted_bam} | grep ^SN | cut -f 2- > ${sample_id}.sorted.bam.stats.tsv    
//     """

// }

process ALIGNMENT_COVERAGES {

    tag "ALIGNMENT_COVERAGES ${sample_id}"
    publishDir "${params.outdir}/alignment-coverages/", pattern: "${sample_id}*.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("*.tsv"), emit: alignment_coverages

    script:
    """
    aligncov.py -i ${sorted_bam} -o ${sample_id}
    """

}

/*
-------------------------------------------------------------------------------
    WORKFLOW
-------------------------------------------------------------------------------
*/

workflow {

    // targets_ch.view()
    // reads_ch.view()

    FASTQC(reads_ch)

    TRIMMOMATIC(reads_ch)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log.collect{it[1]}.ifEmpty([]))
    MULTIQC(ch_multiqc_files.collect())
   
    BWA_INDEX(targets_ch)

    BWA_ALIGN(BWA_INDEX.out.indexed_targets.combine(TRIMMOMATIC.out.trimmed_reads))

    // ALIGNMENT_STATS(BWA_ALIGN.out.aligned_targets)

    ALIGNMENT_COVERAGES(BWA_ALIGN.out.aligned_targets)

}

/*
-------------------------------------------------------------------------------
    EXECUTION SUMMARY
-------------------------------------------------------------------------------
*/

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