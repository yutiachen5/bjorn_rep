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
        final_tsv_translated = Channel.empty()

        ref_id_ch = Channel
            .fromPath(params.ref_file)
            .splitText()
            .filter { it.startsWith(">") }
            .map    { it.replace(">", "").split()[0].trim() }
            .first() // returns the first element in the list, raises error if empty

        tmp_ch = combined_fasta_ch.splitFasta(by: params.chunk_size, file: true, elem: 2)
            .map { f -> tuple(f.baseName, f) }

        alignment_sam = MINIMAP(tmp_ch).alignment_sam
        alignment_fasta = GOFASTA_ALIGNMENT(alignment_sam).alignment_fasta
        

        if (params.gofasta) {
            GOFASTA_VARIANTS(alignment_fasta.combine(ref_id_ch))

            GOFASTA_VARIANTS.out.aa_changes_csv 
                | GOFASTA_CONVERT

            final_tsv = GOFASTA_CONVERT.out.mutations_tsv
                .collectFile(
                    storeDir: "${params.outdir}",
                    name: 'mutations.tsv',
                    newLine: false,
                    skip: 1,
                    keepHeader: true,
                )

            // collect alignemnt and aa_changes files for sanity check
            final_alignment = alignment_fasta
                .map{cid, ali -> ali}
                .collectFile(
                    storeDir: "${params.outdir}",
                    name: 'alignment.fasta',
                    newLine: false,
                    skip: 1*2,
                    keepHeader: true,
                )
            final_aa_changes = GOFASTA_VARIANTS.out.aa_changes_csv
                .map{cid, aa_ch -> aa_ch}
                .collectFile(
                    storeDir: "${params.outdir}",
                    name: 'aa_changes.csv',
                    newLine: false,
                    skip: 1,
                    keepHeader: true,
                )
        } else {
            CALL_MUTATION(alignment_fasta.combine(ref_id_ch))
            mutations_tsv = CALL_MUTATION.out.mutations_tsv
            del_helper_tsv = CALL_MUTATION.out.del_helper_tsv

            final_tsv = CALL_MUTATION.out.mutations_tsv
                .map{ cid, mut -> mut}
                .collectFile(
                    storeDir: "${params.outdir}",
                    name: 'mutations.tsv',
                    newLine: false,
                    skip: 1,
                    keepHeader: true
                )

            // collect helper files for sanity check
            final_helper_tsv = CALL_MUTATION.out.del_helper_tsv
                .map{ cid, del -> del}
                .collectFile(
                    storeDir: "${params.outdir}",
                    name: 'del_helper.tsv',
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
                        final_alignment = alignment_fasta
                            .map{cid, ali -> ali}
                            .collectFile(
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
                        .join(mutations_tsv)
                        .join(del_helper_tsv)
                        .map{ cid, ali, mut, del -> tuple(ali, mut, del)}

                TRANSLATE_MUTATIONS(query_ids_ch, trans_input_ch, ref_id_ch, n_ref)

                final_tsv_translated = TRANSLATE_MUTATIONS.out.mutations_tsv
                    .groupTuple()
                    .flatMap { qid, files ->
                        files.collect { f -> tuple("${qid}_mutations.tsv", f) }
                    }
                    .collectFile(
                        storeDir: "${params.outdir}",
                        name: { it[0] },
                        keepHeader: true,
                        skip: 1,
                        newLine: false
                    )
            }
        }


    emit:
        final_tsv
        final_tsv_translated
}


