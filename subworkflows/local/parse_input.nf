// Modified from https://raw.githubusercontent.com/nf-core/ampliseq/4e48b7100302e2576ac1be2ccc7d464253e9d20e/subworkflows/local/parse_input.nf
workflow PARSE_INPUT {
    take:
    input // folder
    pattern

    main:
    error_message = "\nCannot find any reads matching: \"${input}${pattern}\"\n"
    error_message += "Please revise the input folder (\"--input\"): \"${input}\"\n"
    error_message += "and the input file pattern (\"--pattern\"): \"${pattern}\"\n"
    error_message += "*Please note: Path needs to be enclosed in quotes!*\n"
    error_message += "For more info, please consult the pipeline documentation.\n"
    
    //Get files - paired end
    Channel
        .fromFilePairs( input + pattern, size: 2 )
        .ifEmpty { error("${error_message}") }
        .map { name, reads ->
                def meta = [:]
                meta.id = name.toString().indexOf("_") != -1 ? name.toString().take(name.toString().indexOf("_")) : name
                [ meta, reads ] }
        .set { ch_reads }
    
    //Check whether all sampleID = meta.id are unique
    ch_reads
        .map { meta, reads -> [ meta.id ] }
        .toList()
        .subscribe {
            if( it.size() != it.unique().size() ) {
                ids = it.take(10);
                error("Please review data input, sample IDs are not unique! First IDs are $ids")
            }
        }

    emit:
    reads = ch_reads
}
