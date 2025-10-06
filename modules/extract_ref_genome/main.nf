// extract the ref genome from reference file

process EXTRACT_REFERENCE_GENOME {
    input:
    path ref_genome

    output:
    path "ref.fasta", emit: ref_genome

    script:
    """
    grep -A 1 "^>" ${ref_genome} > ref.fasta
    """
}