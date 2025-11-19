#!/usr/bin/env nextflow

include { EXTRACT_IDS } from '../modules/extract_data/main.nf'
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
        final_tsv_translated = Channel.empty()

        ref_id_ch = Channel
            .fromPath(params.ref_file)
            .splitText()
            .filter { it.startsWith(">") }
            .map    { it.replace(">", "").split()[0].trim() }
            .first() // returns the first element in the list, raises error if empty

        tmp_ch = combined_fasta_ch.splitFasta(by: params.chunk_size, file: true, elem: 2)

        alignment_sam = MINIMAP(tmp_ch).alignment_sam
        alignment_fasta = GOFASTA_ALIGNMENT(alignment_sam).alignment_fasta // the alignment file contains all ref genome
        

        if (params.gofasta) {
            GOFASTA_VARIANTS(alignment_fasta, ref_id_ch)

            GOFASTA_VARIANTS.out.aa_changes_csv 
                | GOFASTA_CONVERT

            final_tsv = GOFASTA_CONVERT.out.mutations_tsv.collectFile(
                storeDir: "${params.outdir}",
                name: 'mutations.tsv',
                newLine: false,
                skip: 1,
                keepHeader: true,
            )

            // collect alignemnt files for sanity check
            final_alignment = alignment_fasta.collectFile(
                storeDir: "${params.outdir}",
                name: 'alignment.fasta',
                newLine: false,
                skip: 1*2,
                keepHeader: true,
            )
        } else {
            mutation_calling_out = CALL_MUTATION(alignment_fasta, ref_id_ch)
            mutations_tsv = mutation_calling_out.mutations_tsv
            del_helper_tsv = mutation_calling_out.del_helper_tsv

            final_tsv = mutations_tsv.collectFile(
                storeDir: "${params.outdir}",
                name: 'mutations.tsv',
                newLine: false,
                skip: 1,
                keepHeader: true
            )

            if (params.translate_mutations) {
                query_ids_ch = Channel
                    .fromPath(params.query_ref_file)
                    .splitText()
                    .filter { it.startsWith(">") }
                    .map    { it.replace(">", "").trim().split()[0] } // multiple ids

                def n_ref = query_ids_ch.toList().size() + ref_id_ch.toList().size()

                // collect alignemnt files for sanity check
                query_ids_ch.toList().then { query_list ->
                    ref_id_ch.toList().then { ref_list ->
                        int n_ref_int = query_list.size() + ref_list.size()   

                        final_alignment = alignment_fasta.collectFile(
                            storeDir: "${params.outdir}",
                            name: 'alignment.fasta',
                            newLine: false,
                            skip: n_ref_int*2,
                            keepHeader: true,
                        )
                    }
                }




                trans_input_ch =
                    alignment_fasta
                        .merge(mutations_tsv)
                        .merge(del_helper_tsv)
                TRANSLATE_MUTATIONS(query_ids_ch, trans_input_ch, ref_id_ch, n_ref)

                // TODO: collect tsv files from all chunks by query id

                // final_tsv = TRANSLATE_MUTATIONS.out.mutations_tsv.collectFile(
                //     storeDir: "${params.outdir}",
                //     name: 'mutations.tsv',
                //     newLine: false,
                //     skip: 1,
                //     keepHeader: true
                // )
            }
        }


    emit:
        final_tsv
        final_tsv_translated
}


