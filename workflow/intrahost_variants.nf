#!/usr/bin/env nextflow

include { SRATOOLS_PREFETCH } from '../modules/nf-core/prefetch/main.nf'
// include { SRATOOLS_FASTERQDUMP } from ''

workflow INTRAHOST_VARIANTS_WORKFLOW {
    take:
    metadata_ch
    reference_ch
    gff_ch
    variant_threshold
    variant_min_depth

    main:

    // prefetch, .map using multiple statements needs curly braces
    prefetch_input = sra_sequences_ch
        .map{
            meta ->
            def id = meta.id
            def sra_id = meta.sra_accession
            [[id: meta.id, sra_accession:sra_id], sra_id]
        }
    SRATOOLS_PREFETCH(prefetch_input)

    // fasteradump
    // SRATOOLS_FASTERQDUMP(SRATOOLS_PREFETCH.out.sra)

    // fastp

    //bwa index

    // bwa mem

    // samtools iview

    //samtools sort streaming ??

}