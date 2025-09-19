nextflow.enable.dsl=2

// variant calling using go_fasta, sam --> csv 
process GOFASTA_VARIANTS {
    input:
    path alignement_sam

    output:
    path alignment_fa, emit: aa_changes

    script:
    """
    gofasta sam variants \\
    -s ${alignment_sam}
    --append-snps \\
    --append-codons \\
    -o aa_changes.csv
    """

}

// file format conversion. csv --> tsv
process GOFASTA_CONVERT {
    input:

    output:

    script:
    """
    
    """

}