#!/usr/bin/env nextflow

include { MINIMAP } from '../modules/alignment/main.nf'
include { MAFFT } from '../modules/alignment/main.nf'
include { GOFASTA_ALIGNMENT } from '../modules/alignment/main.nf'
include { GOFASTA_VARIANTS } from '../modules/mutation_calling/main.nf'
include { GOFASTA_CONVERT } from '../modules/mutation_calling/main.nf'
include { CALL_MUTATION } from '../modules/mutation_calling/main.nf'
include { TRANSLATE_MUTATIONS } from '../modules/mutation_translation/main.nf'

workflow MUTATIONS_ANALYSIS_WORKFLOW {
    take:
        combined_fasta_ch

    main:
        alignment_sam = MINIMAP(combined_fasta_ch).alignment_sam

        // alignment_fasta = MAFFT(combined_fasta_ch).mafft_fasta
        // mutations_tsv = CALL_MUTATION(alignment_fasta).mutations_tsv

        if (params.gofasta) {
            alignment_fasta = GOFASTA_ALIGNMENT(alignment_sam).alignment_fasta
            mutations_csv = GOFASTA_VARIANTS(alignment_fasta).aa_changes_csv
            mutations_tsv = GOFASTA_CONVERT(mutations_csv).mutations_tsv
        } else {
            alignment_fasta = GOFASTA_ALIGNMENT(alignment_sam).alignment_fasta.collect()
            mutations_tsv = CALL_MUTATION(alignment_fasta).mutations_tsv
        }

        if (params.translate_mutations) {
            // Translate mutations
            if (!params.query_id || params.query_id.toString().trim() == '') {
                exit 1, "Query ID of reference genome not specified for translation!"
            } else {
                def query_ids = []

                if (params.query_id instanceof String) {
                    // Handle any separators: comma, space, or both
                    query_ids = params.query_id
                        .replaceAll(/[\[\]"]/, '')     // remove brackets or quotes like ["A","B","C"]
                        .split(/[,\s]+/)               // split by comma or any whitespace
                        .collect { it.trim() }         // trim each
                        .findAll { it }                // remove empties
                }
                else if (params.query_id instanceof List) {
                    // Already a list (e.g. --query_id '[A,B,C]')
                    query_ids = params.query_id.collect { it.toString().trim() }
                }
                else {
                    // Fallback in rare cases
                    query_ids = [ params.query_id.toString().trim() ]
                }
                query_id_ch = query_ids ? Channel.fromList(query_ids) : Channel.empty()
            }

            query_id_ch
                .combine(alignment_fasta)
                .combine(mutation_tsv)
                | TRANSLATE_MUTATIONS
            // TRANSLATE_MUTATIONS(mutations_tsv, alignment_fasta)
        }

    emit:
        mutations_tsv
}


