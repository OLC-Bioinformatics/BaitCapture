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
    IMPORT LOCAL MODULES
-------------------------------------------------------------------------------
*/

// TODO: Transition local modules to nf-core where applicable

// TODO:
// include { ALIGNMENT_COVERAGES as ALIGNMENT_COVERAGES_BWA } from 'modules/local/alignment_coverages'
// include { ALIGNMENT_COVERAGES as ALIGNMENT_COVERAGES_KMA } from 'modules/local/alignment_coverages'
include { ALIGNMENT_COVERAGES_BWA       } from './modules/local/alignment_coverages_bwa'
include { ALIGNMENT_COVERAGES_KMA       } from './modules/local/alignment_coverages_kma'
include { BWA_ALIGN_TARGETS             } from './modules/local/bwa_align_targets'
include { BWA_HOST_DECONTAMINATION      } from './modules/local/bwa_host_decontamination'
// TODO:
// include { BWA_INDEX as BWA_INDEX_HOST } from 'modules/local/bwa_index'
// include { BWA_INDEX as BWA_INDEX_TARGETS } from 'modules/local/bwa_index'
include { BWA_INDEX_HOST                } from './modules/local/bwa_index_host'
include { BWA_INDEX_TARGETS             } from './modules/local/bwa_index_targets'
// TODO:
// include { CONVERT_SAM_TO_SORTED_BAM as CONVERT_SAM_TO_SORTED_BAM_BWA } from 'modules/local/convert_sam_to_sorted_bam'
// include { CONVERT_SAM_TO_SORTED_BAM as CONVERT_SAM_TO_SORTED_BAM_KMA } from 'modules/local/convert_sam_to_sorted_bam'
include { CONVERT_SAM_TO_SORTED_BAM_BWA } from './modules/local/convert_sam_to_sorted_bam_bwa'
include { CONVERT_SAM_TO_SORTED_BAM_KMA } from './modules/local/convert_sam_to_sorted_bam_kma'
include { DECONTAMINATION_STATS         } from './modules/local/decontamination_stats'
include { FASTQC                        } from './modules/local/fastqc'
include { KMA_ALIGN_TARGETS             } from './modules/local/kma_align_targets'
include { KMA_INDEX_TARGETS             } from './modules/local/kma_index_targets'
include { MULTIQC                       } from './modules/local/multiqc'
include { OBTAIN_DECONTAMINATED_READS   } from './modules/local/obtain_decontaminated_reads'
include { TRIMMOMATIC                   } from './modules/local/trimmomatic'

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

        BWA_INDEX_TARGETS(targets_ch)

        KMA_INDEX_TARGETS(targets_ch)

        BWA_INDEX_HOST(host_ch)

        BWA_HOST_DECONTAMINATION(BWA_INDEX_HOST.out.indexed_host.combine(TRIMMOMATIC.out.trimmed_reads))

        OBTAIN_DECONTAMINATED_READS(BWA_HOST_DECONTAMINATION.out.host_aligned_sam)

        DECONTAMINATION_STATS(
            TRIMMOMATIC.out.trimmed_reads.map{
                key, val -> [key, val[0]]
            }.join(OBTAIN_DECONTAMINATED_READS.out.decontaminated_reads.map{
                key, val -> [key, val[0]]
            })
        )

        BWA_ALIGN_TARGETS(BWA_INDEX_TARGETS.out.indexed_targets.combine(OBTAIN_DECONTAMINATED_READS.out.decontaminated_reads))

        KMA_ALIGN_TARGETS(KMA_INDEX_TARGETS.out.indexed_targets.combine(OBTAIN_DECONTAMINATED_READS.out.decontaminated_reads))

        CONVERT_SAM_TO_SORTED_BAM_BWA(BWA_ALIGN_TARGETS.out.targets_aligned_sam)

        CONVERT_SAM_TO_SORTED_BAM_KMA(KMA_ALIGN_TARGETS.out.targets_aligned_sam)

        ALIGNMENT_COVERAGES_BWA(CONVERT_SAM_TO_SORTED_BAM_BWA.out.aligned_sorted_bam)

        ALIGNMENT_COVERAGES_KMA(CONVERT_SAM_TO_SORTED_BAM_KMA.out.aligned_sorted_bam)

    }
    else {

        BWA_INDEX_TARGETS(targets_ch)

        BWA_ALIGN_TARGETS(BWA_INDEX_TARGETS.out.indexed_targets.combine(TRIMMOMATIC.out.trimmed_reads))

        KMA_INDEX_TARGETS(targets_ch)

        KMA_ALIGN_TARGETS(KMA_INDEX_TARGETS.out.indexed_targets.combine(TRIMMOMATIC.out.trimmed_reads))

        CONVERT_SAM_TO_SORTED_BAM_BWA(BWA_ALIGN_TARGETS.out.targets_aligned_sam)

        CONVERT_SAM_TO_SORTED_BAM_KMA(KMA_ALIGN_TARGETS.out.targets_aligned_sam)

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