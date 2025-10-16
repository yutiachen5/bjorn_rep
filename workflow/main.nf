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
        if (params.sampling) {
            fasta_files = SAMPLING(params.fasta_dir).selected_fasta
        } else {
            fasta_files = GET_ALL_FASTA(params.fasta_dir).all_fasta
        }
        
        log.info "Extracting mutations..."
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

        if (params.query_id == null) {
            query_id_ch = Channel.empty()
        } else {
            Channel
                .fromList(params.query_id)
                .set { query_id_ch }
        }

        if (params.translate_mutations) {
            QUERY_GENOMES(mutations_tsv, query_id_ch)
        }

    // emit:
    // mutations_tsv
}