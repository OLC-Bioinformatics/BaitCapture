#!/usr/bin/env nextflow

/*
-------------------------------------------------------------------------------
    BaitCapture
-------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include {validateParameters; paramsHelp; paramsSummaryLog} from 'plugin/nf-validation'

// Define parameters
params.reads = "*R{1,2}.fastq.gz"
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

if (params.host) {
    host_ch = Channel.fromPath("${params.host}", checkIfExists: true)
} else {
    host_ch = []
}

/*
-------------------------------------------------------------------------------
    PROCESSES
-------------------------------------------------------------------------------
*/

process FASTQC {

    conda "bioconda::fastqc=0.11.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'quay.io/biocontainers/fastqc:0.11.9--0' }"

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
        'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2' }"

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

    conda "bioconda::multiqc=1.15"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.15--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0' }"

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

    conda "bioconda::bwa=0.7.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
        'quay.io/biocontainers/bwa:0.7.17--hed695b0_7' }"

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

    conda "bioconda::bwa=0.7.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--he4a0461_11' :
        'quay.io/biocontainers/bwa:0.7.17--he4a0461_11' }"

    cpus = 8
    tag "BWA_ALIGN ${sample_id}"

    input:
    tuple path(targets), path("*"), val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bwa.aligned.sam"), emit: targets_aligned_sam

    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem -t ${task.cpus} \$INDEX ${reads} > ${sample_id}.bwa.aligned.sam
    """

}

process KMA_INDEX {

    conda "bioconda::kma=1.4.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.9--he4a0461_2' :
        'quay.io/biocontainers/kma:1.4.9--he4a0461_2' }"

    tag "KMA_INDEX ${targets}"

    input:
    path targets

    output:
    tuple path(targets), path("*"), emit: indexed_targets

    script:
    """
    kma index -i ${targets} -o targets.kma
    """

}

process KMA_ALIGN {

    conda "bioconda::kma=1.4.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kma:1.4.9--he4a0461_2' :
        'quay.io/biocontainers/kma:1.4.9--he4a0461_2' }"

    cpus = 8
    tag "KMA_ALIGN ${sample_id}"   

    input:
    tuple path(indexed_targets), path("*"), val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.kma.aligned.sam"), emit: targets_aligned_sam

    script:
    """
    kma -mem_mode -ex_mode -1t1 -vcf \
        -t ${task.cpus} \
        -ipe ${reads} \
        -t_db targets.kma \
        -o ${sample_id}.kma.aligned.temp \
        -sam \
        > ${sample_id}.kma.aligned.sam
    """

}

process CONVERT_SAM_TO_SORTED_BAM_BWA {

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_1' :
        'quay.io/biocontainers/samtools:1.17--hd87286a_1' }"

    cpus = 8
    tag "CONVERT_SAM_TO_SORTED_BAM_BWA ${sample_id}"
    publishDir "${params.outdir}/bwa/", pattern: "${sample_id}.bwa.aligned.sorted.bam", mode: 'copy'

    input:
    tuple val(sample_id), path(targets_aligned_sam)

    output:
    tuple val(sample_id), path("${sample_id}.bwa.aligned.sorted.bam"), emit: aligned_sorted_bam

    script:
    """
    samtools view --threads ${task.cpus} -S -b ${sample_id}.bwa.aligned.sam > ${sample_id}.bwa.aligned.bam
    samtools sort --threads ${task.cpus} ${sample_id}.bwa.aligned.bam > ${sample_id}.bwa.aligned.sorted.bam
    """

}

process CONVERT_SAM_TO_SORTED_BAM_KMA {

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_1' :
        'quay.io/biocontainers/samtools:1.17--hd87286a_1' }"

    cpus = 8
    tag "CONVERT_SAM_TO_SORTED_BAM_KMA ${sample_id}"
    publishDir "${params.outdir}/kma/", pattern: "${sample_id}.kma.aligned.sorted.bam", mode: 'copy'

    input:
    tuple val(sample_id), path(targets_aligned_sam)

    output:
    tuple val(sample_id), path("${sample_id}.kma.aligned.sorted.bam"), emit: aligned_sorted_bam

    script:
    """
    samtools view --threads ${task.cpus} -S -b ${sample_id}.kma.aligned.sam > ${sample_id}.kma.aligned.bam
    samtools sort --threads ${task.cpus} ${sample_id}.kma.aligned.bam > ${sample_id}.kma.aligned.sorted.bam
    """

}

process ALIGNMENT_COVERAGES_BWA {

    conda "bioconda::aligncov=0.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aligncov:0.0.2--pyh7cba7a3_0' :
        'quay.io/biocontainers/aligncov:0.0.2--pyh7cba7a3_0' }"

    tag "ALIGNMENT_COVERAGES_BWA ${sample_id}"
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

process ALIGNMENT_COVERAGES_KMA {

    conda "bioconda::aligncov=0.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aligncov:0.0.2--pyh7cba7a3_0' :
        'quay.io/biocontainers/aligncov:0.0.2--pyh7cba7a3_0' }"

    tag "ALIGNMENT_COVERAGES_KMA ${sample_id}"
    publishDir "${params.outdir}/kma/", pattern: "${sample_id}*.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("*.tsv"), emit: alignment_coverages

    script:
    """
    aligncov -i ${sorted_bam} -o ${sample_id}.kma
    """

}

process INDEX_HOST {

    conda "bioconda::bwa=0.7.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--hed695b0_7' :
        'quay.io/biocontainers/bwa:0.7.17--hed695b0_7' }"

    tag "INDEX_HOST ${host}"

    input:
    path host

    output:
    tuple path(host), path("*"), emit: indexed_host

    script:
    """
    bwa index ${host}
    """

}

process HOST_DECONTAMINATION {

    conda "bioconda::bwa=0.7.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.17--he4a0461_11' :
        'quay.io/biocontainers/bwa:0.7.17--he4a0461_11' }"

    cpus = 8
    tag "HOST_DECONTAMINATION ${sample_id}"

    input:
    tuple path(host), path("*"), val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.host.aligned.sam"), emit: host_aligned_sam

    script:
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem -t ${task.cpus} \$INDEX ${reads} > ${sample_id}.host.aligned.sam
    """

}

process OBTAIN_DECONTAMINATED_READS {

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_1' :
        'quay.io/biocontainers/samtools:1.17--hd87286a_1' }"

    cpus = 8
    tag "OBTAIN_DECONTAMINATED_READS ${sample_id}"
    publishDir "${params.outdir}/decontaminated-reads/", pattern: "${sample_id}*.fastq", mode: 'copy'

    input:
    tuple val(sample_id), path(host_aligned_sam)

    output:
    tuple val(sample_id), path("${sample_id}.host.aligned.sorted.bam"), emit: host_aligned_sorted_bam
    tuple val(sample_id), path("${sample_id}_R{1,2}.decontam.fastq"), emit: decontaminated_reads

    script:
    """
    # Convert SAM to sorted BAM
    samtools view --threads ${task.cpus} -b -f 0x0004 -f 0x0008 -f 0x0001 ${sample_id}.host.aligned.sam | samtools sort --threads ${task.cpus} > ${sample_id}.host.aligned.sorted.bam
    # Convert sorted BAM to FASTQ
    samtools fastq --threads ${task.cpus} -1 ${sample_id}_R1.decontam.fastq -2 ${sample_id}_R2.decontam.fastq -s shouldnotexist.fq ${sample_id}.host.aligned.sorted.bam
    """

}

process DECONTAMINATION_STATS {

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--hd87286a_1' :
        'quay.io/biocontainers/samtools:1.17--hd87286a_1' }"

    tag "DECONTAMINATION_STATS ${sample_id}"
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
   
   // If a host genome is provided in CLI arguments:
    if (params.host) {

        BWA_INDEX(targets_ch)

        KMA_INDEX(targets_ch)

        INDEX_HOST(host_ch)

        HOST_DECONTAMINATION(INDEX_HOST.out.indexed_host.combine(TRIMMOMATIC.out.trimmed_reads))

        OBTAIN_DECONTAMINATED_READS(HOST_DECONTAMINATION.out.host_aligned_sam)

        DECONTAMINATION_STATS(
            TRIMMOMATIC.out.trimmed_reads.map{
                key, val -> [key, val[0]]
            }.join(OBTAIN_DECONTAMINATED_READS.out.decontaminated_reads.map{
                key, val -> [key, val[0]]
            })
        )

        BWA_ALIGN(BWA_INDEX.out.indexed_targets.combine(OBTAIN_DECONTAMINATED_READS.out.decontaminated_reads))

        KMA_ALIGN(KMA_INDEX.out.indexed_targets.combine(OBTAIN_DECONTAMINATED_READS.out.decontaminated_reads))

        CONVERT_SAM_TO_SORTED_BAM_BWA(BWA_ALIGN.out.targets_aligned_sam)

        CONVERT_SAM_TO_SORTED_BAM_KMA(KMA_ALIGN.out.targets_aligned_sam)

        ALIGNMENT_COVERAGES_BWA(CONVERT_SAM_TO_SORTED_BAM_BWA.out.aligned_sorted_bam)

        ALIGNMENT_COVERAGES_KMA(CONVERT_SAM_TO_SORTED_BAM_KMA.out.aligned_sorted_bam)

    }
    else {

        BWA_INDEX(targets_ch)

        BWA_ALIGN(BWA_INDEX.out.indexed_targets.combine(TRIMMOMATIC.out.trimmed_reads))

        KMA_INDEX(targets_ch)

        KMA_ALIGN(KMA_INDEX.out.indexed_targets.combine(TRIMMOMATIC.out.trimmed_reads))

        CONVERT_SAM_TO_SORTED_BAM_BWA(BWA_ALIGN.out.targets_aligned_sam)

        CONVERT_SAM_TO_SORTED_BAM_KMA(KMA_ALIGN.out.targets_aligned_sam)

        ALIGNMENT_COVERAGES_BWA(CONVERT_SAM_TO_SORTED_BAM_BWA.out.aligned_sorted_bam)

        ALIGNMENT_COVERAGES_KMA(CONVERT_SAM_TO_SORTED_BAM_KMA.out.aligned_sorted_bam)

    } 

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