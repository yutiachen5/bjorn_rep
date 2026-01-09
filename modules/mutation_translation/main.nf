#!/usr/bin/env nextflow

process TRANSLATE_MUTATIONS { 
    tag "${chunk_id}"
    label "process_high"

    input:
    each query_id
    tuple val(chunk_id), path(alignment_fasta), path(mutations_tsv), path(del_helper_tsv)
    val(ref_id)
    val(n_ref)

    output:
    tuple val(chunk_id), val(query_id), path("${query_id}_mutations.tsv"), emit: mutations_tsv
    path("versions.yml"), emit: versions

    script:
    """
    translate_mutations.py \
        -a ${alignment_fasta} \
        -m ${mutations_tsv} \
        -d ${del_helper_tsv} \
        --gff ${params.gff_file} \
        --query ${query_id} \
        --bg ${ref_id} \
        --n_ref ${n_ref} \
        --region ${query_id} \
        -o ${query_id}_mutations.tsv 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${query_id}_mutations.tsv
    touch versions.yml
    """
}
