// alignment using minimap2

process MINIMAP {
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

process GOFASTA_ALIGNMENT {
    input:
    path alignment_sam
    path ref_genome
    
    output:
    path "alignment.fasta", emit: alignment_fasta

    script:
    """
    gofasta sam toMultiAlign -s ${alignment_sam} \\
                             -o alignment_samples.fasta

    cat ${ref_genome} alignment_samples.fasta > alignment_ref_samples.fasta 

    seqkit seq -w 0 alignment_ref_samples.fasta > alignment.fasta                      
    """
}


