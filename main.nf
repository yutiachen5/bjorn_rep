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
}

workflow {
    log.info "Start consensus genome alignment.."

    validateParameters()

    // segments_ch = Channel.empty()

    GENERAL_WORKFLOW()
}