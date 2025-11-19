// randomly sample 100 samples from the consensus genome or get all fasta files under the input dir

process SAMPLING {

    output:
    path "combined.fasta", emit: combined_fasta

    script:
    """
    mkdir -p ${params.outdir}

    sampling.py --fasta_dir ${params.fasta_dir} \
                --nsamples ${params.nsamples} \
                -s ${params.seed} \
                -o fasta_files.txt

    xargs cat < fasta_files.txt > combined.fasta

    cp -p combined.fasta ${params.outdir}/combined.fasta
    """
    
}

process GET_ALL_FASTA {

    output:
    path "combined.fasta", emit: combined_fasta

    script:
    """
    mkdir -p ${params.outdir}

    find ${params.fasta_dir} -type f \\( -name "*.fasta" -o -name "*.fa" \\) -exec realpath {} \\; > fasta_files.txt

    xargs cat < fasta_files.txt > combined.fasta

    cp -p combined.fasta ${params.outdir}/combined.fasta
    """

}

process EXTRACT_IDS {

    input:
        path ref_file
        path query_file

    output:
        stdout into ids_ch

    script:
    """
    mapfile -t ref_ids < <(grep '^>' ${ref_file} | sed 's/^>//' | awk '{print \$1}')
    n_bg_ref=\${#ref_ids[@]}

    if [[ \$n_bg_ref -ne 1 ]]; then
        echo "ERROR: Reference file '${ref_file}' contains \$n_bg_ref sequences; exactly 1 required." >&2
        exit 1
    fi

    ref_id="\${ref_ids[0]}"

    mapfile -t query_ids_arr < <(grep '^>' ${query_file} | sed 's/^>//' | awk '{print \$1}')
    n_query_ref=\${#query_ids_arr[@]}

    query_ids="\${query_ids_arr[0]}"
    for id in "\${query_ids_arr[@]:1}"; do
        query_ids="\${query_ids},\${id}"
    done

    n_ref=\$(( 1 + n_query_ref ))

    printf "%s\n" "$ref_id" "$query_ids" "$n_ref"
    """
}
