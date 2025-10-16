#!/usr/bin/env nextflow

include { MAFFT } from '../modules/alignment/main.nf'
include { MINIMAP } from '../modules/alignment/main.nf'
include { TRANSLATE_MUTATIONS } from '../modules/mutation_translation/main.nf'
include { GOFASTA_ALIGNMENT } from '../modules/alignment/main.nf'

workflow QUERY_GENOMES {
    take:
        mutation_ch
        query_id_ch

    main:
        // mafft_ch = MAFFT(params.ref_file, params.query).mafft_fasta.collect()

        alignment_sam = MINIMAP(params.query).alignment_sam
        alignment_fasta = GOFASTA_ALIGNMENT(alignment_sam).alignment_fasta.collect()

        TRANSLATE_MUTATIONS(alignment_fasta, mutation_ch.collect(), query_id_ch)

    // emit:

}