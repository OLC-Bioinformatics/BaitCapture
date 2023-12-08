include { BWAMEM2_INDEX as BWAMEM2_BUILD } from '../../modules/nf-core/bwamem2/index/main'
include { BWAMEM2_MEM as BWAMEM2_ALIGN   } from '../../modules/nf-core/bwamem2/mem/main'

workflow BWAMEM2_ALIGN_READS {

    take:
    ch_targets
    ch_final_reads

    main:
    ch_versions = Channel.empty()

    // Build index for gene target database
    BWAMEM2_BUILD(ch_targets)
    ch_indexed_targets = BWAMEM2_BUILD.out.index
    ch_versions = ch_versions.mix(BWAMEM2_BUILD.out.versions.first())

    // Align final reads to indexed target database
    BWAMEM2_ALIGN(ch_final_reads, ch_indexed_targets.collect(), true)
    ch_sorted_bam = BWAMEM2_ALIGN.out.bam
    ch_versions = ch_versions.mix(BWAMEM2_ALIGN.out.versions.first())

    emit:
    sorted_bam = ch_sorted_bam
    versions = ch_versions

}