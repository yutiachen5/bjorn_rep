process MINIMAP {
    input:
    path query_seq

    output:
    path "alignment.sam", emit: alignment_sam

    script:
    """
    cat ${params.ref_file} ${query_seq} > query.fasta
    minimap2 -a \\
             -x asm20 \\
             --score-N=0 \\
             -N 0 \\
             -p 0.8 \\
             ${params.ref_file} query.fasta > alignment.sam
    """
}

process MAFFT {
    input:
    path ref_genome
    path query_genome

    output:
    path "mafft.fasta", emit: mafft_fasta

    script:
    """
    cat ${ref_genome} ${query_genome} > cat_tmp.fasta
    
    mafft cat_tmp.fasta > mafft_alignment.fasta

    seqkit seq -w 0 mafft_alignment.fasta > mafft.fasta
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
                             -o alignment.fasta                
    """
}


