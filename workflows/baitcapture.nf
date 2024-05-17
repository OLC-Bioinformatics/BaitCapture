/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { samplesheetToList; paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-schema'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowBaitcapture.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ALIGNCOV                      } from '../modules/local/aligncov/main'
include { BWAMEM2_ALIGN_READS           } from '../subworkflows/local/bwamem2_align_reads'
include { MERGE_MAPPING_RESULTS         } from '../modules/local/merge_mapping_results'
include { PARSE_INPUT                   } from '../subworkflows/local/parse_input'
include { BWA_ALIGN_READS               } from '../subworkflows/local/bwa_align_reads'
include { KMA_ALIGN_READS               } from '../subworkflows/local/kma_align_reads'
include { SUMMARIZE_STATS               } from '../modules/local/summarize_stats'
include { VALIDATE_TARGET_METADATA      } from '../modules/local/validate_target_metadata'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BAM_STATS_SAMTOOLS as BAM_STATS_SAMTOOLS_HOST          } from '../subworkflows/nf-core/bam_stats_samtools/main'
include { BAM_STATS_SAMTOOLS as BAM_STATS_SAMTOOLS_TARGETS       } from '../subworkflows/nf-core/bam_stats_samtools/main'
include { BWAMEM2_HOST_REMOVAL_MEM as BWAMEM2_HOST_REMOVAL_ALIGN } from '../modules/local/bwamem2_host_mem/main'
include { BWAMEM2_INDEX as BWAMEM2_HOST_REMOVAL_BUILD            } from '../modules/nf-core/bwamem2/index/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                            } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTQC as FASTQC_RAW                                   } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_PREPROCESSED                          } from '../modules/nf-core/fastqc/main'
include { FASTQSCAN as FASTQSCAN_RAW                             } from '../modules/nf-core/fastqscan/main'
include { FASTQSCAN as FASTQSCAN_PREPROCESSED                    } from '../modules/nf-core/fastqscan/main'
include { MOSDEPTH                                               } from '../modules/nf-core/mosdepth/main'
include { MULTIQC                                                } from '../modules/nf-core/multiqc/main'
include { FASTP                                                  } from '../modules/nf-core/fastp/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_HOST                  } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TARGETS               } from '../modules/nf-core/samtools/index/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow BAITCAPTURE {

    ch_versions = Channel.empty()

    if (params.input) {
            // Argument is the name of the parameter which specifies the samplesheet, i.e. params.input = 'input'
            // [[id:SRR14739083], [s3://sra-pub-src-5/SRR14739083/HI.4968.007.UDI0015_i7---UDI0015_i5.D9_SP3-0_mock_R1.fastq.gz.1, s3://sra-pub-src-5/SRR14739083/HI.4968.007.UDI0015_i7---UDI0015_i5.D9_SP3-0_mock_R2.fastq.gz.1]]
            ch_reads = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json")).map{
                meta, fastq_1, fastq_2 ->
                def reads = [fastq_1, fastq_2]
                return [meta, reads]
            }
    } else if (params.input_folder) {
        PARSE_INPUT(params.input_folder, params.extension)
        ch_reads = PARSE_INPUT.out.reads
    } else {
        error("One of `--input` or `--input_folder` must be provided!")
    }    

    if (params.targets) {
        ch_targets = Channel.fromPath(params.targets, checkIfExists: true).map{
            def meta = [:]
            meta.id = 'targets'
            [meta, it]
        }
    } else {
        error("A FASTA file `--targets` must be provided!")
    }

    if (params.host) {
        ch_host = Channel.fromPath(params.host, checkIfExists: true).map{
            def meta = [:]
            meta.id = 'host'
            [meta, it]
        }
    }

    if (params.adapters) {
        ch_adapters = Channel.fromPath(params.adapters, checkIfExists: true)
    }

    if (params.target_metadata) {
        ch_target_metadata = Channel.fromPath(params.target_metadata, checkIfExists: true)
        VALIDATE_TARGET_METADATA(ch_targets, ch_target_metadata)
        ch_versions = ch_versions.mix(VALIDATE_TARGET_METADATA.out.versions.first())
    }

    //
    // MODULE: FASTQC_RAW
    //
    FASTQC_RAW(ch_reads)
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    //
    // MODULE: FASTP
    //
    if (!params.skip_trimming) {
        if (params.adapters) {
            FASTP(ch_reads, ch_adapters.collect(), [], [])
            ch_trimmed_reads = FASTP.out.reads
        } else {
            FASTP(ch_reads, [], [], [])
            ch_trimmed_reads = FASTP.out.reads
        }
        ch_fastp_log = FASTP.out.log
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    }

    //
    // MODULE: BWAMEM2_HOST_REMOVAL_BUILD
    //
    if (params.host) {
        BWAMEM2_HOST_REMOVAL_BUILD(ch_host)
        ch_host_index = BWAMEM2_HOST_REMOVAL_BUILD.out.index
        ch_versions = ch_versions.mix(BWAMEM2_HOST_REMOVAL_BUILD.out.versions.first())
    }

    //
    // MODULE: BWAMEM2_HOST_REMOVAL_ALIGN
    //

    if (params.host) {
        if (!params.skip_trimming) {
            // Trimming + host decontamination
            BWAMEM2_HOST_REMOVAL_ALIGN(ch_trimmed_reads, ch_host_index.collect())
            ch_preprocessed_reads = BWAMEM2_HOST_REMOVAL_ALIGN.out.reads
            ch_sorted_bam_host = BWAMEM2_HOST_REMOVAL_ALIGN.out.bam
        } else {
            // Host decontamination only
            BWAMEM2_HOST_REMOVAL_ALIGN(ch_reads, ch_host_index.collect())
            ch_preprocessed_reads = BWAMEM2_HOST_REMOVAL_ALIGN.out.reads
            ch_sorted_bam_host = BWAMEM2_HOST_REMOVAL_ALIGN.out.bam
        }
    } else if (!params.skip_trimming) {
            // Trimming only
            ch_preprocessed_reads = ch_trimmed_reads
    }

    //
    // MODULE: SAMTOOLS_INDEX_HOST
    //
    if (params.host) {
        SAMTOOLS_INDEX_HOST(ch_sorted_bam_host)
        ch_sorted_bam_bai_host = SAMTOOLS_INDEX_HOST.out.bai
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_HOST.out.versions.first())
    }

    //
    // SUBWORKFLOW: BAM_STATS_SAMTOOLS_TARGETS
    //
    if (params.host){
        BAM_STATS_SAMTOOLS_HOST(ch_sorted_bam_host.join(ch_sorted_bam_bai_host, by: [0]), ch_targets.collect())
        ch_idxstats = BAM_STATS_SAMTOOLS_HOST.out.idxstats
        ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS_HOST.out.versions.first())
    }

    //
    // MODULE: FASTQC_PREPROCESSED
    //
    if (params.host || !params.skip_trimming) {
        FASTQC_PREPROCESSED(ch_preprocessed_reads)
        ch_versions = ch_versions.mix(FASTQC_PREPROCESSED.out.versions.first())
    }

    /*
    ========================================================================
    Obtain input and preprocessed read counts
    ========================================================================
    */

    //
    // MODULE: FASTQSCAN_RAW
    //
    FASTQSCAN_RAW(ch_reads)
    ch_fastqscan_raw = FASTQSCAN_RAW.out.json
    ch_versions = ch_versions.mix(FASTQSCAN_RAW.out.versions.first())

    //
    // MODULE: FASTQSCAN_PREPROCESSED
    //
    if (params.host) {
        FASTQSCAN_PREPROCESSED(ch_preprocessed_reads)
        ch_fastqscan_preprocessed = FASTQSCAN_PREPROCESSED.out.json
        ch_versions = ch_versions.mix(FASTQSCAN_RAW.out.versions.first())
    }

    /*
    ========================================================================
    Align (preprocessed) reads to targets
    ========================================================================
    */

   if (params.host || !params.skip_trimming) {
        ch_final_reads = ch_preprocessed_reads
    } else {
        // No trimming or host decontamination
        ch_final_reads = ch_reads
    }

    if (params.aligner == 'bwamem2') {
        BWAMEM2_ALIGN_READS(ch_targets, ch_final_reads)
        ch_sorted_bam_targets = BWAMEM2_ALIGN_READS.out.sorted_bam
        ch_versions = ch_versions.mix(BWAMEM2_ALIGN_READS.out.versions)
    } else if (params.aligner == 'kma') {
        KMA_ALIGN_READS(ch_targets, ch_final_reads)
        ch_sorted_bam_targets = KMA_ALIGN_READS.out.sorted_bam
        ch_kma_res = KMA_ALIGN_READS.out.res
        ch_versions = ch_versions.mix(KMA_ALIGN_READS.out.versions)
    } else if (params.aligner == 'bwa') {
        BWA_ALIGN_READS(ch_targets, ch_final_reads)
        ch_sorted_bam_targets = BWA_ALIGN_READS.out.sorted_bam
        ch_versions = ch_versions.mix(BWA_ALIGN_READS.out.versions)
    }

    //
    // MODULE: ALIGNCOV
    //
    ALIGNCOV(ch_sorted_bam_targets)
    ch_aligncov_stats = ALIGNCOV.out.stats

    //
    // MODULE: SAMTOOLS_INDEX_TARGETS
    //
    SAMTOOLS_INDEX_TARGETS(ch_sorted_bam_targets)
    ch_sorted_bam_bai_targets = SAMTOOLS_INDEX_TARGETS.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_TARGETS.out.versions.first())

    //
    // SUBWORKFLOW: BAM_STATS_SAMTOOLS_TARGETS
    //
    BAM_STATS_SAMTOOLS_TARGETS(ch_sorted_bam_targets.join(ch_sorted_bam_bai_targets, by: [0]), ch_targets.collect())
    ch_idxstats = BAM_STATS_SAMTOOLS_TARGETS.out.idxstats
    ch_stats = BAM_STATS_SAMTOOLS_TARGETS.out.stats
    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS_TARGETS.out.versions.first())

    /*
    ========================================================================
    Analysis of mapping results
    ========================================================================
    */

    //
    // MODULE: MOSDEPTH
    //
    ch_mosdepth_in = ch_sorted_bam_targets.join(ch_sorted_bam_bai_targets, by: [0]).map{ meta, bam, bai -> [meta, bam, bai, []] }
    MOSDEPTH(ch_mosdepth_in, ch_targets.collect())
    ch_mosdepth_summary = MOSDEPTH.out.summary_txt
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    //
    // MODULE: MERGE_MAPPING_RESULTS
    //
    if (params.aligner == 'kma') {
        ch_merge_mapping_results_in = ch_aligncov_stats
        .combine(ch_idxstats, by: [0])
        .combine(ch_kma_res, by: [0])
    } else {
        ch_merge_mapping_results_in = ch_aligncov_stats
        .combine(ch_idxstats, by: [0])
        .combine(ch_aligncov_stats.map{ meta, kma_res -> [meta, []] }, by: [0])
    }
    if (params.target_metadata) {
        MERGE_MAPPING_RESULTS(ch_merge_mapping_results_in, ch_target_metadata.collect())
    } else {
        MERGE_MAPPING_RESULTS(ch_merge_mapping_results_in, [])
    }
    ch_versions = ch_versions.mix(MERGE_MAPPING_RESULTS.out.versions.first())

    //
    // MODULE: SUMMARIZE_STATS
    //
    if (!params.skip_trimming && !params.host) {
        ch_summarize_stats_in = ch_fastqscan_raw
        .combine(ch_fastqscan_raw.map{ meta, fastqscan_preprocessed -> [meta, []] }, by: [0])
        .combine(ch_stats, by: [0])
        .combine(ch_fastp_log, by: [0])
    } else if (params.skip_trimming && !params.host) {
        ch_summarize_stats_in = ch_fastqscan_raw
        .combine(ch_fastqscan_raw.map{ meta, fastqscan_preprocessed -> [meta, []] }, by: [0])
        .combine(ch_stats, by: [0])
        .combine(ch_fastqscan_raw.map{ meta, fastp_log -> [meta, []] }, by: [0])
    } else if (!params.skip_trimming && params.host) {
        ch_summarize_stats_in = ch_fastqscan_raw
        .combine(ch_fastqscan_preprocessed, by: [0])
        .combine(ch_stats, by: [0])
        .combine(ch_fastp_log, by: [0])
    } else if (params.skip_trimming && params.host) {
        ch_summarize_stats_in = ch_fastqscan_raw
        .combine(ch_fastqscan_preprocessed, by: [0])
        .combine(ch_stats, by: [0])
        .combine(ch_fastqscan_raw.map{ meta, fastp_log -> [meta, []] }, by: [0])
    }
    SUMMARIZE_STATS(ch_summarize_stats_in)
    ch_versions = ch_versions.mix(SUMMARIZE_STATS.out.versions.first())

    /*
    ========================================================================
    Workflow cleanup
    ========================================================================
    */

    //
    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MULTIQC
    //
    workflow_summary    = WorkflowBaitcapture.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowBaitcapture.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]}.ifEmpty([]))
    if (params.host || !params.skip_trimming) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PREPROCESSED.out.zip.collect{it[1]}.ifEmpty([]))
    }
    if (!params.skip_trimming) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    }
    // ch_multiqc_files = ch_multiqc_files.mix(BAM_STATS_SAMTOOLS_TARGETS.out.idxstats.collect{it[1]}.ifEmpty([]))


    // TODO: Add process outputs for MultiQC input here

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
