#!/usr/bin/env nextflow

include { MINIMAP } from '../modules/alignment/main.nf'
include { GOFASTA_ALIGNMENT } from '../modules/alignment/main.nf'
include { GOFASTA_VARIANTS } from '../modules/mutation_calling/main.nf'
include { GOFASTA_CONVERT } from '../modules/mutation_calling/main.nf'
include { CALL_MUTATION } from '../modules/mutation_calling/main.nf'

workflow MUTATIONS_ANALYSIS_WORKFLOW {
    take:
        combined_fasta_ch

    main:
        alignment_sam = MINIMAP(params.ref_file, combined_fasta_ch).alignment_sam

        mutations_csv = GOFASTA_VARIANTS(alignment_sam).aa_changes_csv
        mutations_tsv = GOFASTA_CONVERT(mutations_csv).mutations_tsv

        // alignment_fasta = GOFASTA_ALIGNMENT(alignment_sam, params.ref_file).alignment_fasta.collect()
        // mutations_tsv = CALL_MUTATION(alignment_fasta).mutations_tsv

    emit:
        mutations_tsv
}


