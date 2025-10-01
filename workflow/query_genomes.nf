#!/usr/bin/env nextflow

include { MAFFT } from '../modules/alignment/main.nf'
include { TRANSLATE_MUTATIONS } from '../modules/mutation_translation/main.nf'

workflow QUERY_GENOMES {
    take:
        mutation_ch
        query_id_ch

    main:
        // alignment file between target and query genome
        mafft_ch = MAFFT(params.ref_file, params.query).mafft_fasta .collect()
        mutation_ch = mutation_ch .collect()
        TRANSLATE_MUTATIONS(mafft_ch, mutation_ch, query_id_ch)

    // emit:


}