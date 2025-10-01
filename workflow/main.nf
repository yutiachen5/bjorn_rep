#!/usr/bin/env nextflow

include { SAMPLING } from '../modules/sampling/main.nf'
include { VARIANT_ANALYSIS_WORKFLOW } from './variant_analysis.nf'
include { INTRAHOST_VARIANTS_WORKFLOW } from './intrahost_variants.nf' 

workflow EXTRACT_MUTATIONS {
    take:
    metadata_ch

    main:
        SAMPLING(params.fasta_dir)
        selected_fasta = SAMPLING.out.selected_fasta

        if (params.run_intrahost_variants) { // false currently
            log.info "Running intrahost variant analysis"

            sra_sequences_ch = metadata_ch
                .splitCsv(header: true, sep:"\t")
                .filter { row -> row.SraAccessions && row.SraAccessions != ""}
                .map { row -> 
                    def meta = [id: row.Accession, accession: row.Accession, sra_accession: row.SraAccessions]
                    return meta 
                }

            // Intrahost variant
            INTRAHOST_VARIANTS(
                sra_sequences_ch,
                Channel.value(params.ref_file),
                Channel.fromPath(params.ivar_gff),  
                params.variant_threshold,
                params.variant_min_depth
            )
        } else {
            log.info "Skipping intrahost variant analysis"
        }

        // EXTRACT_REFERENCE_GENOME(params.ref_file)
        // ref_genome = EXTRACT_REFERENCE_GENOME.out.ref_file_genome

        log.info "Running variant analysis using consensus sequence on the whole genome"
        fasta_ch = selected_fasta
                    .splitText()
                    .map{ filename -> 
                        def trimed_filename = filename.trim()
                        file(trimed_filename).normalize() 
                        } 
        combined_fasta_ch = fasta_ch.collectFile() { filename ->
            ["combined.fasta", file(filename).text]
        }

        mutations_tsv = VARIANT_ANALYSIS_WORKFLOW(combined_fasta_ch).mutations_tsv

    emit:
    mutations_tsv
}