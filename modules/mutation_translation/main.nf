#!/usr/bin/env nextflow

process TRANSLATE_MUTATIONS {    
    input:
    path alignment_fasta
    path mutations_tsv
    val ref_id
    val query_id

    output:
    path "${query_id}_mutations.tsv", emit: mutations_tsv

    script:
    """
    ref_id=\$(grep '^>' ${params.ref_file} | head -n 1 | sed 's/^>//; s/ .*//')

    translate_mutations.py \
        -a ${alignment_fasta} \
        -m ${mutations_tsv} \
        --gff ${params.gff_file} \
        --query ${query_id} \
        --bg ${ref_id} \
        --region ${params.region} \
        -o mutations.tsv 
        
    cp -p ${query_id}_mutations.tsv ${params.outdir}/${query_id}_mutations.tsv
    """
}