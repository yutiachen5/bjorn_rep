#!/usr/bin/env nextflow

process TRANSLATE_MUTATIONS {
    input:
    path mafft_fasta
    path mutations_tsv
    val query_id

    output:
    path "${query_id}_mutations.tsv", emit: mutations_tsv

    script:
    """
    translate_mutations.py \
        -a ${mafft_fasta} \
        -o mutations.tsv \
        -m ${mutations_tsv} \
        --query ${query_id} \
        --ref ${params.ref_id}
    cp -p ${query_id}_mutations.tsv ${params.outdir}/${query_id}_mutations.tsv
    """
}