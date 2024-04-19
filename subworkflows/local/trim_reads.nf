include { FASTP       } from '../../modules/nf-core/fastp/main'
include { TRIMMOMATIC } from '../../modules/nf-core/trimmomatic/main'

workflow TRIM_READS {

    take:
    ch_reads

    main:
    ch_versions = Channel.empty()

    // Use FASTP by default for trimming unless --use_trimmomatic provided
    if (params.use_trimmomatic) {
        TRIMMOMATIC(ch_reads)
        ch_trimmed_reads = TRIMMOMATIC.out.trimmed_reads
        ch_json = Channel.empty()
        ch_versions = ch_versions.mix(TRIMMOMATIC.out.versions.first())
    } else {
        FASTP(ch_reads, [], [], [])
        ch_trimmed_reads = FASTP.out.reads
        ch_json = FASTP.out.json
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    }

    emit:
    json          = ch_json
    trimmed_reads = ch_trimmed_reads
    versions      = ch_versions

}