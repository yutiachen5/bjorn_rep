#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// main file to take the input

include { GENERAL_WORKFLOW } from './workflow/main.nf'

def validateParameters() {
    if (!params.ref) {
        exit 1, "Reference genome file not specified!"
    }
    if (!params.fasta_path) {
        exit 1, "Input FASTA file not specified!"
    }
    if (!params.gb_dir) {
        exit 1, "Input GB directory not specified!"
    }
    // if (!params.metadata) {
    //     exit 1, "Metadata not specified!"
    // }
}

workflow {
    log.info "Start consensus genome alignment.."

    validateParameters()

    // segments_ch = Channel.empty()

    metadata_ch = Channel.fromPath(params.metadata)

    GENERAL_WORKFLOW(metadata_ch)
}