#!/usr/bin/env nextflow

include { SAMPLING } from '../modules/extract_data/main.nf'
include { GET_ALL_FASTA } from '../modules/extract_data/main.nf'
include { MUTATIONS_ANALYSIS_WORKFLOW } from './mutations_analysis.nf'
include { INTRAHOST_VARIANTS_WORKFLOW } from './intrahost_variants.nf' 
include { QUERY_GENOMES } from './query_genomes.nf'


workflow EXTRACT_MUTATIONS {
    take:
    metadata_ch

    main:
        // Extract mutations
        if (params.sampling) {
            fasta_files = SAMPLING(params.fasta_dir).selected_fasta
        } else {
            fasta_files = GET_ALL_FASTA(params.fasta_dir).all_fasta
        }
        
        fasta_ch = fasta_files
                    .splitText()
                    .map{ filename -> 
                        def trimed_filename = filename.trim()
                        file(trimed_filename).normalize() 
                        } 
        combined_fasta_ch = fasta_ch.collectFile() { filename ->
            ["combined.fasta", file(filename).text]
        }

        mutations_tsv = MUTATIONS_ANALYSIS_WORKFLOW(combined_fasta_ch).mutations_tsv

        // Translate mutations
        if (!params.query_id || params.query_id.toString().trim() == '') {
            // No input provided
            query_id_ch = Channel.empty()
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


        if (params.translate_mutations) {
            QUERY_GENOMES(mutations_tsv, query_id_ch)
        }

    // emit:
    // mutations_tsv
}