// alignment using minimap2

process DO_MINIMAP {
    input:
    path ref_genome
    path query_seq

    output:
    path "alignment.sam", emit: alignment_sam

    script:
    """
    minimap2 -a \\
             -x asm20 \\
             --score-N=0 \\
             -N 0 \\
             -p 0.8 \\
             ${ref_genome} ${query_seq} > alignment.sam
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


