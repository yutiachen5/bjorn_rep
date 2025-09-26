// extract the ref genome from reference file

process EXTRACT_REFERENCE_GENOME {
    input:
    path ref_genome_path

    output:
    path "ref.fasta", emit: ref_genome

    script:
    """
    echo ${ref_genome_path}
    grep -A 1 "^>" ${ref_genome_path} > ref.fasta
    """
}