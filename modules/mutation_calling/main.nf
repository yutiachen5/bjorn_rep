// mutation calling using go_fasta, sam --> csv 
process GOFASTA_SAM_VARIANTS {
    input:
    path alignment_sam

    output:
    path "aa_changes.csv", emit: aa_changes_csv

    script:
    """
    mkdir -p ${params.outdir}

    gofasta sam variants -s ${alignment_sam} \\
                         -a ${params.gff_file} \\
                         --append-snps \\
                         -o aa_changes.csv

    cp -p aa_changes.csv ${params.outdir}/aa_changes.csv
    """
}

// mutation calling using go_fasta, fasta --> csv 
process GOFASTA_VARIANTS {
    input:
    path alignment_fasta
    val ref_id

    output:
    path "aa_changes.csv", emit: aa_changes_csv

    script:
    """
    mkdir -p ${params.outdir}

    gofasta variants --msa ${alignment_fasta} \\
                     --reference ${ref_id} \\
                     -a ${params.gff_file} \\
                     --append-snps \\
                     -o aa_changes.csv
                     
    cp -p aa_changes.csv ${params.outdir}/aa_changes.csv
    """
}

// file format conversion, csv --> tsv
process GOFASTA_CONVERT {
    input:
    path aa_changes_csv

    output:
    path "mutations.tsv", emit: mutations_tsv

    script:
    """
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
    val ref_id
    
    output:
    path "mutations.tsv", emit: mutations_tsv
    path "del_helper.tsv", emit: del_helper_tsv

    script:
    """
    mkdir -p ${params.outdir}

    cat ${params.ref_file} ${params.query_ref_file} > ref_all.fasta

    mutation_calling.py -a ${alignment_fasta} \\
                        --gff ${params.gff_file} \\
                        --region ${params.region} \\
                        --ref_file ref_all.fasta \\
                        --ref_id ${ref_id} \\
                        -o mutations.tsv

    cp -p mutations.tsv ${params.outdir}/mutations.tsv
    cp -p del_helper.tsv ${params.outdir}/del_helper.tsv
    """
}