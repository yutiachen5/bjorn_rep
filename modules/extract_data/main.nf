// randomly sample 100 samples from the consensus genome or get all fasta files under the input dir

process SAMPLING {
    input:
    val fasta_dir

    output:
    path "fasta_files.txt", emit: selected_fasta

    script:
    if (params.l) {
        """
        sampling.py --fasta_dir ${params.fasta_dir} \
                    -l ${params.l} \
                    --lineage_file ${params.lineage_file} \
        """
    } else {
        """
        sampling.py --fasta_dir ${params.fasta_dir} \
                    -n 100 \
        """
    }
}

process GET_ALL_FASTA {
    input:
    val fasta_dir

    output:
    path "fasta_files.txt", emit: all_fasta

    script:
    """
    find ${fasta_dir} -type f \\( -name "*.fasta" -o -name "*.fa" \\) -exec realpath {} \\; > fasta_files.txt
    """

}