
process MINIMAP {
    input:
    path query_seq

    output:
    path "alignment.sam", emit: alignment_sam

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

    cp -p all_seq.fasta ${params.outdir}/all_seq.fasta
    cp -p alignment.sam ${params.outdir}/alignment.sam
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
    input:
    path alignment_sam
    
    output:
    path "alignment.fasta", emit: alignment_fasta

    script:
    """
    gofasta sam toMultiAlign -s ${alignment_sam} \\
                             -r ${params.ref_file} \\
                             -o alignment.fasta 
    cp -p alignment.fasta ${params.outdir}/alignment.fasta               
    """
}


