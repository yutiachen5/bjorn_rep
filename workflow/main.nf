#!/usr/bin/env nextflow

include { SAMPLING } from '../modules/extract_data/main.nf'
include { GET_ALL_FASTA } from '../modules/extract_data/main.nf'
include { MUTATIONS_ANALYSIS_WORKFLOW } from './mutations_analysis.nf'

workflow EXTRACT_MUTATIONS {

    main:
        if (params.sampling) {
            combined_fasta_ch = SAMPLING().combined_fasta
        } else {
            combined_fasta_ch = GET_ALL_FASTA().combined_fasta
        }

        mutations_tsv = MUTATIONS_ANALYSIS_WORKFLOW(combined_fasta_ch)
}
