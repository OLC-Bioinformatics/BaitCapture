include { BWA_INDEX as BWA_BUILD } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM as BWA_ALIGN   } from '../../modules/nf-core/bwa/mem/main'

workflow BWA_ALIGN_READS {

    take:
    ch_targets
    ch_final_reads

    main:
    ch_versions = Channel.empty()

    // Build index for gene target database
    BWA_BUILD(ch_targets)
    ch_indexed_targets = BWA_BUILD.out.index
    ch_versions = ch_versions.mix(BWA_BUILD.out.versions.first())

    // Align final reads to indexed target database
    BWA_ALIGN(ch_final_reads, ch_indexed_targets.collect(), [[id:'no_fasta'], []], true)
    ch_sorted_bam = BWA_ALIGN.out.bam
    ch_versions = ch_versions.mix(BWA_ALIGN.out.versions.first())

    emit:
    sorted_bam = ch_sorted_bam
    versions = ch_versions

}