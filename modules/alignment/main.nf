nextflow.enable.dsl=2

// alignment using minimap2

process DO_MINIMAP {
    input:
    path ref_genome
    path combined_fasta

    output:
    path "alignment.sam", emit: alignment_sam

    script:
    """
    ls -l ${ref_genome}
    ls -l ${combined_fasta}
    minimap2 -a \\
             -x asm20 \\
             --score-N=0 \\
             -N 0 \\
             -p 0.8 \\
             ${ref_genome} ${combined_fasta} > alignment.sam
    """

}
