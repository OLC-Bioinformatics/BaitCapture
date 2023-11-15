/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { fromSamplesheet; paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

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

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { PARSE_INPUT                   } from '../subworkflows/local/parse_input'
include { ALIGNCOV as ALIGNCOV_BWAMEM2  } from '../modules/local/aligncov/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                        } from '../modules/nf-core/fastqc/main'
include { MULTIQC                       } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { TRIMMOMATIC                   } from '../modules/nf-core/trimmomatic/main'
include { BWAMEM2_INDEX as BWAMEM2_BUILD } from '../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_ALIGN  } from '../modules/nf-core/bwamem2/mem/main'
include { BOWTIE2_REMOVAL_BUILD as BOWTIE2_HOST_REMOVAL_BUILD } from '../modules/local/bowtie2_removal_build'
include { BOWTIE2_REMOVAL_ALIGN as BOWTIE2_HOST_REMOVAL_ALIGN } from '../modules/local/bowtie2_removal_align'

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
            ch_reads = Channel.fromSamplesheet('input').map{
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
        ch_host = Channel.fromPath(params.host, checkIfExists: true)
    }

    //
    // MODULE: FASTQC
    //
    FASTQC(ch_reads)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: TRIMMOMATIC
    //
    if (!params.skip_trimmomatic) {
        TRIMMOMATIC(ch_reads)
        ch_trimmed_reads = TRIMMOMATIC.out.trimmed_reads
        ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())
    }

    //
    // MODULE: BOWTIE2_HOST_REMOVAL_BUILD
    //
    if (params.host) {
        BOWTIE2_HOST_REMOVAL_BUILD(ch_host)
        ch_host_index = BOWTIE2_HOST_REMOVAL_BUILD.out.index
        ch_versions = ch_versions.mix(BOWTIE2_HOST_REMOVAL_BUILD.out.versions.first())
    }

    //
    // MODULE: BOWTIE2_HOST_REMOVAL_ALIGN
    //
    if (params.host) {
        if (!params.skip_trimmomatic) {
            BOWTIE2_HOST_REMOVAL_ALIGN(ch_trimmed_reads, ch_host_index)
            ch_trimmed_reads_decontaminated = BOWTIE2_HOST_REMOVAL_ALIGN.out.reads
        } else {
            BOWTIE2_HOST_REMOVAL_ALIGN(ch_reads, ch_host_index)
            ch_reads_decontaminated = BOWTIE2_HOST_REMOVAL_ALIGN.out.reads
        }
        ch_versions = ch_versions.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.versions.first())
    }

    //
    // MODULE: BWAMEM2_BUILD
    //
    BWAMEM2_BUILD(ch_targets)
    ch_indexed_targets = BWAMEM2_BUILD.out.index
    ch_versions = ch_versions.mix(BWAMEM2_BUILD.out.versions.first())

    //
    // MODULE: BWAMEM2_ALIGN
    //
    if (!params.skip_trimmomatic) {

        if (params.host) {
            BWAMEM2_ALIGN(ch_trimmed_reads_decontaminated, ch_indexed_targets, 'sort')
            ch_bwa_sorted_bam = BWAMEM2_ALIGN.out.bam
        } else {
            BWAMEM2_ALIGN(ch_trimmed_reads, ch_indexed_targets, 'sort')
            ch_bwa_sorted_bam = BWAMEM2_ALIGN.out.bam            
        }

    } else {

        if (params.host) {
            BWAMEM2_ALIGN(ch_reads_decontaminated, ch_indexed_targets, 'sort')
            ch_bwa_sorted_bam = BWAMEM2_ALIGN.out.bam
        } else {
            BWAMEM2_ALIGN(ch_treads, ch_indexed_targets, 'sort')
            ch_bwa_sorted_bam = BWAMEM2_ALIGN.out.bam
        }

    }
    ch_versions = ch_versions.mix(BWAMEM2_ALIGN.out.versions.first())

    //
    // MODULE: ALIGNCOV_BWAMEM2
    //
    ALIGNCOV_BWAMEM2(ch_bwa_sorted_bam)

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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.log.collect{it[1]}.ifEmpty([]))
    if (!params.skip_trimmomatic) {
        ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.log.collect{it[1]}.ifEmpty([]))
    }
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
