// mutation calling using go_fasta, sam --> csv 
process GOFASTA_SAM_VARIANTS {
    input:
    path alignment_sam

    output:
    path "aa_changes.csv", emit: aa_changes_csv

    script:
    """
    gofasta sam variants -s ${alignment_sam} \\
                         -a ${params.gff_file} \\
                         --append-snps \\
                         --append-codons \\
                         -o aa_changes.csv
    """
}

// mutation calling using go_fasta, fasta --> csv 
process GOFASTA_VARIANTS {
    tag "${chunk_id}"
    label "process_medium"

    input:
    tuple val(chunk_id), path(alignment_fasta), val(ref_id)

    output:
    tuple val(chunk_id), path("aa_changes.csv"), emit: aa_changes_csv
    path("versions.yml"), emit: versions

    script:
    """
    gofasta variants --msa ${alignment_fasta} \\
                     --reference ${ref_id} \\
                     -a ${params.gff_file} \\
                     --append-snps \\
                     --append-codons \\
                     -o aa_changes.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gofasta: \$(gofasta --version 2>&1 | sed 's/gofasta //g')
    END_VERSIONS
    """

    stub:
    """
    touch aa_changes.csv
    touch versions.yml
    """
}

// file format conversion, csv --> tsv
process GOFASTA_CONVERT {
    tag "${chunk_id}"
    label "process_low"

    input:
    tuple val(chunk_id), path(aa_changes_csv)

    output:
    path("mutations.tsv"), emit: mutations_tsv
    path("versions.yml"), emit: versions

    script:
    """
    gofasta_converter.py --csv_path ${aa_changes_csv} \\
                        --reference ${params.ref_file} \\
                        --region ${params.region} \\
                        --output mutations.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch mutations.tsv
    touch versions.yml
    """
}

// manual mutation calling, include missing values  and X as invalid aa
process CALL_MUTATION {
    tag "${chunk_id}"
    label 'process_high'

    input:
    tuple val(chunk_id), path(alignment_fasta), val(ref_id)
    
    output:
    tuple val(chunk_id), path("${ref_id}_mutations.tsv"), emit: mutations_tsv
    tuple val(chunk_id), path("del_helper.tsv"), emit: del_helper_tsv
    path("versions.yml"), emit: versions

    script:
    """
    cat ${params.ref_file} ${params.query_ref_file} > ref_all.fasta

    mutation_calling.py -a ${alignment_fasta} \\
                        --gff ${params.gff_file} \\
                        --region ${params.region} \\
                        --ref_file ref_all.fasta \\
                        --ref_id ${ref_id} \\
                        -o ${ref_id}_mutations.tsv 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${ref_id}_mutations.tsv
    touch del_helper.tsv
    touch versions.yml
    """
}

// stack mutation files
process CONCAT_MUTATIONS {
    label "process_low"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path mut_files, stageAs: 'in_*'

    output:
    path 'mutations.tsv'

    script:
    """
    set -euo pipefail
    shopt -s nullglob

    files=(in_*)
    head -n 1 "\${files[0]}" > mutations.tsv
    for f in "\${files[@]}"; do
      tail -n +2 "\$f" >> mutations.tsv
    done
    """
}
