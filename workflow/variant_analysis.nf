#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DO_MINIMAP } from '../modules/alignment/main.nf'
include { GOFASTA_VARIANTS } from '../modules/variant_calling/main.nf'
include { GOFASTA_CONVERT } from '../modules/variant_calling/main.nf'

workflow VARIANT_ANALYSIS_WORKFLOW {
    take:
        ref
        selected_fasta

    main:
        DO_MINIMAP(ref, selected_fasta)
        res_sam = DO_MINIMAP.out.alignment_sam

        GOFASTA_VARIANTS(res_sam)
        res_csv = GOFASTA_VARIANTS.out.aa_changes_csv

        GOFASTA_CONVERT(res_csv)
        res_tsv = GOFASTA_CONVERT.out.mutations_tsv

}


