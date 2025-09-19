nextflow.enable.dsl=2

include { SAMPLING } from '../modules/sampling/main.nf'
include { EXTRACT_REFERENCE_GENOME } from '../modules/extract_ref_genome/main.nf'

include { VARIANT_ANALYSIS_WORKFLOW } from './variant_analysis.nf'

workflow GENERAL_WORKFLOW {
    main:
        SAMPLING(params.fasta_path)
        selected_fasta = SAMPLING.out.selected_fasta

        EXTRACT_REFERENCE_GENOME(params.ref)
        ref_genome = EXTRACT_REFERENCE_GENOME.out.ref_genome

        if (params.run_intrahost_variants) { // false currently
            log.info "Running intrahost variant analysis"
        } else {
            log.info "Skipping intrahost variant analysis"
        }

        log.info "Running variant analysis using consensus sequence on the whole genome"
        fasta_ch = selected_fasta
                    .splitText()
                    .map{ filename -> file("${params.fasta_path}/${filename}") } // make each string as a single file
        fasta_ch.view()

        VARIANT_ANALYSIS_WORKFLOW(ref, selected_segments)

        res_sam = DO_MINIMAP(ref, selected_segments)
        res_csv = GO_FASTA(res_sam)
}