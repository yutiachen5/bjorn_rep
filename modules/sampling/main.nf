nextflow.enable.dsl=2

// randomly sample 100 segments from the consensus genome

process SAMPLING {
    input:
    val consensus_genome_path

    output:
    path "selected_segments.txt", emit: selected_fasta

    script:
    """
    # have to limit to first 1000 lines due to memory restriction (exitcode 143), may find other ways later
    ls ${consensus_genome_path} | head -n 1000 | shuf -n 100 > selected_segments.txt
    """
}