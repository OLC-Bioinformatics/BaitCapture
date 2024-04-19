include { FASTP       } from '../../modules/nf-core/fastp/main'
include { TRIMMOMATIC } from '../../modules/nf-core/trimmomatic/main'

workflow TRIM_READS {

    take:
    ch_reads

    main:
    ch_versions = Channel.empty()

    TRIMMOMATIC(ch_reads)
    ch_trimmed_reads = TRIMMOMATIC.out.trimmed_reads
    ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())

    emit:
    trimmed_reads = ch_trimmed_reads
    versions = ch_versions

}