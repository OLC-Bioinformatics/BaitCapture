include { FASTP       } from '../../modules/nf-core/fastp/main'

workflow TRIM_READS {

    take:
    ch_reads

    main:
    ch_versions = Channel.empty()

    FASTP(ch_reads, [], [], [])
    ch_trimmed_reads = FASTP.out.reads
    ch_json = FASTP.out.json
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    emit:
    json          = ch_json
    trimmed_reads = ch_trimmed_reads
    versions      = ch_versions

}