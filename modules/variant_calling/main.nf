nextflow.enable.dsl=2

// variant calling using go_fasta, sam --> csv 
process GOFASTA_VARIANTS {
    input:
    path alignment_sam

    output:
    path "aa_changes.csv", emit: aa_changes_csv

    script:
    """
    gofasta sam variants \\
    -s ${alignment_sam} \\
    -a ${params.gb_dir} \\
    --append-snps \\
    -o aa_changes.csv
    """

}

// file format conversion. csv --> tsv
process GOFASTA_CONVERT {
    input:
    path aa_changes_csv

    output:
    path "mutations.tsv", emit: mutations_tsv

    script:
    """
    echo ${params.outdir}
    gofasta_converter.py --csv_path ${aa_changes_csv} --reference ${params.ref} --region None
    """

}