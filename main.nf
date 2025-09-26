#!/usr/bin/env nextflow

include { GENERAL_WORKFLOW } from './workflow/main.nf'

def validateParameters() {
    if (!params.ref_file) {
        exit 1, "Reference genome file not specified!"
    }
    if (!params.gff_file) {
        exit 1, "Input GFF directory not specified!"
    }
    if (!params.fasta_dir) {
        exit 1, "Input FASTA directory not specified!"
    }
    // if (!params.metadata) {
    //     exit 1, "Metadata not specified!"
    // }
}

workflow {
    log.info "Start consensus genome alignment.."

    validateParameters()

    // segments_ch = Channel.empty()

    // metadata_ch = Channel.fromPath(params.metadata)
    metadata_ch = Channel.empty()

    GENERAL_WORKFLOW(metadata_ch)
}