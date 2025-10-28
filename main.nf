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
    if (params.translate_mutations && !params.query_file) {
        exit 1, "Query genome not specified!"
    }
}

def helpMessage() {
    log.info """
    PB2 DMS Nextflow Pipeline
    -------------------------
    Required parameters:
      --fasta_dir       Path to directory containing sample FASTA files
      --ref_file        Reference FASTA file
      --gff_file        GFF annotation file
      --query           Query FASTA file
      --region          Region name (e.g. PP755596.1)
      --outdir          Output directory

    Optional parameters:
      --ref_id          Reference ID (e.g. 'lcl|PP755596.1_cds_XAJ25426.1_1')
      --query_id        One or more Query IDs, comma- or space-separated.
                        Examples:
                            --query_id A
                            --query_id A,B,C
                            --query_id 'lcl|CY018884.1_cds_ABM21959.1_1,lcl|CY018885.1_cds_ABM21960.1_1'

      --nsamples        Number of random samples (default: 1000)
      -c nf.config      Use custom configuration file
    """
}

if (params.help) {
    helpMessage()
    System.exit(0)
}

workflow {
    log.info "Start consensus genome alignment.."

    validateParameters()

    // segments_ch = Channel.empty()

    // metadata_ch = Channel.fromPath(params.metadata)
    metadata_ch = Channel.empty()

    // EXTRACT_MUTATIONS(metadata_ch)

    mutation_ch = EXTRACT_MUTATIONS(metadata_ch).mutations_tsv

}