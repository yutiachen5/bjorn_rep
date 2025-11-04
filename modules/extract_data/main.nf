// randomly sample 100 samples from the consensus genome or get all fasta files under the input dir

process SAMPLING {

    output:
    path "combined.fasta", emit: combined_fasta

    script:
    """
    mkdir -p ${params.outdir}

    sampling.py --fasta_dir ${params.fasta_dir} \
                --nsamples ${params.nsamples} \
                -o fasta_files.txt

    xargs cat < fasta_files.txt > combined.fasta

    cp -p combined.fasta ${params.outdir}/combined.fasta
    """
    
}

process GET_ALL_FASTA {

    output:
    path "combined.fasta", emit: combined_fasta

    script:
    """
    find ${params.fasta_dir} -type f \\( -name "*.fasta" -o -name "*.fa" \\) -exec realpath {} \\; > fasta_files.txt

    xargs cat < fasta_files.txt > combined.fasta

    cp -p combined.fasta ${params.outdir}/combined.fasta
    """

}