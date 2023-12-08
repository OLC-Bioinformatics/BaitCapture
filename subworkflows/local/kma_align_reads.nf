include { KMA_INDEX                              } from '../../modules/local/kma_index'
include { KMA_ALIGN                              } from '../../modules/local/kma_align'
include { SAMTOOLS_SORT as KMA_SAM_TO_SORTED_BAM } from '../../modules/nf-core/samtools/sort/main'

workflow KMA_ALIGN_READS {

    take:
    ch_targets
    ch_final_reads

    main:
    ch_versions = Channel.empty()

    // Build index for gene target database
    KMA_INDEX(ch_targets)
    ch_indexed_targets = KMA_INDEX.out.index
    ch_versions = ch_versions.mix(KMA_INDEX.out.versions.first())

    // Align final reads to indexed target database
    KMA_ALIGN(ch_final_reads, ch_indexed_targets.collect())
    ch_sam = KMA_ALIGN.out.sam
    ch_versions = ch_versions.mix(KMA_ALIGN.out.versions.first())

    // Convert SAM to sorted BAM
    KMA_SAM_TO_SORTED_BAM(ch_sam)
    ch_sorted_bam = KMA_SAM_TO_SORTED_BAM.out.bam
    ch_versions = ch_versions.mix(KMA_SAM_TO_SORTED_BAM.out.versions.first())

    emit:
    sorted_bam = ch_sorted_bam
    versions = ch_versions

}