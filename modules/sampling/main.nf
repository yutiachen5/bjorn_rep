// randomly sample 100 samples from the consensus genome

process SAMPLING {
    input:
    val fasta_dir

    output:
    path "selected_fasta.txt", emit: selected_fasta

    script:
    if (params.l) {
        """
        sampling.py --fasta_dir ${params.fasta_dir} \
                    -l ${params.l} \
                    --lineage_file ${params.lineage_file}
        """
    } else {
        """
        sampling.py --fasta_dir ${params.fasta_dir} -n 100
        """
    }
}