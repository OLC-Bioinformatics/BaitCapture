include { KMA_INDEX } from '../../modules/local/kma/index'
include { KMA_ALIGN } from '../../modules/local/kma/align'

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
    KMA_ALIGN(ch_final_reads, ch_indexed_targets.collect(), true)
    ch_sorted_bam = KMA_ALIGN.out.bam
    ch_res = KMA_ALIGN.out.res
    ch_versions = ch_versions.mix(KMA_ALIGN.out.versions.first())

    emit:
    sorted_bam = ch_sorted_bam
    res = ch_res
    versions = ch_versions

}