#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { EXTRACT_MUTATIONS } from './workflow/main.nf'

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
    if (params.translate_mutations && !params.query) {
        exit 1, "Query genome not specified!"
    }
}

workflow {
    log.info "Start consensus genome alignment.."

    validateParameters()

    // segments_ch = Channel.empty()

    // metadata_ch = Channel.fromPath(params.metadata)
    if (params.query_id == null) {
        query_id_ch = Channel.empty()
    } else {
        Channel
            .fromList(params.query_id)
            .set { query_id_ch }
    }

    metadata_ch = Channel.empty()

    // EXTRACT_MUTATIONS(metadata_ch)

    mutation_ch = EXTRACT_MUTATIONS(metadata_ch).mutations_tsv

}