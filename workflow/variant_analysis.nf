nextflow.enable.dsl=2

include { DO_MINIMAP, GO_FASTA } from '../modules/alignment/main.nf'
// include { GO_FASTA } from '../modules/alignment/main.nf'

workflow VARIANT_ANALYSIS {
    take:
        ref, segments

    main:
        res_sam = DO_MINIMAP(ref, selected_segments)
        res_csv = GO_FASTA(res_sam)
}


