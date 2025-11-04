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
        query_id_ch
            .combine(alignment_fasta)
            .combine(mutation_ch)
            | TRANSLATE_MUTATIONS

}