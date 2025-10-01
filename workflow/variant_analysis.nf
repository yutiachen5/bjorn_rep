#!/usr/bin/env nextflow

include { DO_MINIMAP } from '../modules/alignment/main.nf'
include { GOFASTA_VARIANTS } from '../modules/variant_calling/main.nf'
include { GOFASTA_CONVERT } from '../modules/variant_calling/main.nf'

workflow VARIANT_ANALYSIS_WORKFLOW {
    take:
        selected_fasta

    main:
        DO_MINIMAP(params.ref_file, selected_fasta)
        res_sam = DO_MINIMAP.out.alignment_sam

        GOFASTA_VARIANTS(res_sam)
        mutations_csv = GOFASTA_VARIANTS.out.aa_changes_csv

        GOFASTA_CONVERT(mutations_csv)
        mutations_tsv = GOFASTA_CONVERT.out.mutations_tsv

    emit:
        mutations_tsv
}


