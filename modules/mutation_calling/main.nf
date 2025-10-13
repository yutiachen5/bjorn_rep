// variant calling using go_fasta, sam --> csv 
process GOFASTA_VARIANTS {
    input:
    path alignment_sam

    output:
    path "aa_changes.csv", emit: aa_changes_csv

    script:
    """
    gofasta sam variants -s ${alignment_sam} \\
                         -a ${params.gff_file} \\
                         --append-snps \\
                         -o aa_changes.csv

    cp -p aa_changes.csv ${params.outdir}/aa_changes.csv
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
    mkdir -p ${params.outdir}

    gofasta_converter.py --csv_path ${aa_changes_csv} \\
                        --reference ${params.ref_file} \\
                        --region ${params.region} \\
                        --output mutations.tsv

    cp -p mutations.tsv ${params.outdir}/mutations.tsv
    """
}

// manual mutation calling, include missing values 
process CALL_MUTATION {
    input:
    path alignment_fasta
    
    output:
    path "mutations.tsv", emit: mutations_tsv

    script:
    """
    mkdir -p ${params.outdir}

    mutation_calling.py -a ${alignment_fasta} \\
                        --gff ${params.gff_file} \\
                        --region ${params.region} \\
                        -o mutations.tsv

    cp -p mutations.tsv ${params.outdir}/mutations.tsv
    """
}