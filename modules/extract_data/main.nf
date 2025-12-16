process SAMPLING {
    label "process_low"

    output:
    path("combined.fasta"), emit: combined_fasta
    path("versions.yml"), emit: versions

    script:
    """
    mkdir -p ${params.outdir}

    sampling.py --fasta_dir ${params.fasta_dir} \
                --nsamples ${params.nsamples} \
                -s ${params.seed} \
                -o fasta_files.txt

    xargs cat < fasta_files.txt > combined.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    echo ">stub_sequence" > combined.fasta
    echo "ACGTACGTACGTACGT" >> combined.fasta
    touch versions.yml
    """
}

process GET_ALL_FASTA {
    label "process_low"

    output:
    path "combined.fasta", emit: combined_fasta

    script:
    """
    mkdir -p ${params.outdir}

    find ${params.fasta_dir} -type f \\( -name "*.fasta" -o -name "*.fa" \\) -exec realpath {} \\; > fasta_files.txt

    xargs cat < fasta_files.txt > combined.fasta
    """

    stub:
    """
    echo ">stub_sequence" > combined.fasta
    echo "ACGTACGTACGTACGT" >> combined.fasta
    """
}
