process MINIMAP {
    input:
    path query_seq

    output:
    path "alignment.sam", emit: alignment_sam

    script:
    """
    mkdir -p ${params.outdir}

    cat ${params.ref_file} ${query_seq} > query.fasta

    minimap2 -a \\
             -x map-ont \\
             --score-N=0 \\
             --sam-hit-only \\
             --secondary=no \\
             ${params.ref_file} query.fasta > alignment.sam

    cp -p query.fasta ${params.outdir}/query.fasta
    cp -p alignment.sam ${params.outdir}/alignment.sam
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
    mkdir -p ${params.outdir}

    gofasta sam toMultiAlign -s ${alignment_sam} \\
                             -r ${params.ref_file} \\
                             -o alignment.fasta 
    cp -p alignment.fasta ${params.outdir}/alignment.fasta               
    """
}


