#!/usr/bin/env nextflow

include { PREFETCH } from ''

workflow INTRAHOST_VARIANTS_WORKFLOWN {
    take:
    metadata_ch,
    reference_ch,
    gff_ch,
    variant_threshold,
    variant_min_depth

    main:

    // prefetch
    prefetch_input = sra_sequences_ch
        .map(
            meta ->
            def id = meta.id
            def sra_id = meta.sra_accession
            [[id: meta.id, sra_accession:sra_id], sra_id]
        )
    SRATOOLS_PREFETCH(prefetch_input)

    // fasteradump

    // fastp

    //bwa index

    // bwa mem

    // samtools iview

    //samtools sort streaming ??

}