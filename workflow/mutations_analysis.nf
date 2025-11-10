#!/usr/bin/env nextflow

include { MINIMAP } from '../modules/alignment/main.nf'
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
        alignment_fasta = GOFASTA_ALIGNMENT(alignment_sam).alignment_fasta // the alignment file contains all ref genome

        ref_id = Channel
                    .fromPath(params.ref_file)
                    .splitText()
                    .map { line ->
                        if (line.startsWith(">")) line.replace(">", "").trim().split()[0]
                    }
                    .filter { it }    
                    .first()

        if (params.gofasta) {
            mutations_csv = GOFASTA_VARIANTS(alignment_fasta, ref_id).aa_changes_csv
            mutations_tsv = GOFASTA_CONVERT(mutations_csv).mutations_tsv
        } else {
            (mutations_tsv, del_helper_tsv) = CALL_MUTATION(alignment_fasta, ref_id)
        }

        if (params.translate_mutations) {
            query_id_ch = Channel
                            .fromPath(params.query_ref_file)
                            .splitText()
                            .map { line ->
                                if (line.startsWith(">")) line.replace(">", "").trim().split()[0]
                            }
                            .filter { it }
            unique_count_ch = query_id_ch.unique().count()+1 // total number of reference geome

            TRANSLATE_MUTATIONS(alignment_fasta, mutations_tsv, del_helper_tsv, ref_id, query_id_ch, unique_count_ch)
        }

    emit:
        mutations_tsv
}


