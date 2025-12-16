process MINIMAP {
    tag "${chunk_id}"
    label "process_medium"

    input:
    tuple val(chunk_id), path(query_seq)

    output:
    tuple val(chunk_id), path("alignment.sam"), emit: alignment_sam
    path("versions.yml"), emit: versions

    script:
    """
    if [ -f "${params.query_ref_file}" ]; then
        cat ${params.ref_file} ${params.query_ref_file} ${query_seq} > all_seq.fasta
    else
        cat ${params.ref_file} ${query_seq} > all_seq.fasta
    fi

    minimap2 -a \\
             -x asm20 \\
             --score-N=0 \\
             --sam-hit-only \\
             --secondary=no \\
             ${params.ref_file} all_seq.fasta > alignment.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1 | sed 's/minimap2 //g')
    END_VERSIONS
    """

    stub:
    """
    touch alignment.sam
    touch versions.yml
    """
}

process MAFFT {
    input:
    path query_seq

    output:
    path "mafft.fasta", emit: mafft_fasta

    script:
    """
    mafft --auto --keeplength --addfragments ${query_seq} ${params.ref_file} | seqkit seq -w 0 > mafft.fasta
    """
}

process GOFASTA_ALIGNMENT {
    tag "${chunk_id}"
    label "process_medium"

    input:
    tuple val(chunk_id), path(alignment_sam)
    
    output:
    tuple val(chunk_id), path("alignment.fasta"), emit: alignment_fasta
    path("versions.yml"), emit: versions

    script:
    """
    gofasta sam toMultiAlign -s ${alignment_sam} \\
                             -r ${params.ref_file} \\
                             -o alignment.fasta 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gofasta: \$(gofasta --version 2>&1 | sed 's/gofasta //g')
    END_VERSIONS
    """

    stub:
    """
    touch alignment.fasta
    touch versions.yml
    """
}


