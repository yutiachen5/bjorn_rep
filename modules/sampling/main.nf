// randomly sample 100 samples from the consensus genome

process SAMPLING {
    input:
    val fasta_dir

    output:
    path "selected_fasta.txt", emit: selected_fasta

    script:
    """
    # echo ${fasta_dir}
    # have to limit to first 1000 lines due to memory restriction (exitcode 143), may find other ways later
    # ls ${fasta_dir} | head -n 1000 | shuf -n 100 > selected_fasta.txt

    sampling.py --fasta_dir ${params.fasta_dir} 
    """
}