nextflow.enable.dsl=2

// alignment using minimap2

process DO_MINIMAP {
    input:
    tuple path(ref_genome), path(fasta_file)

    output:
    path "aa_changes.csv", emit: aa_changes

    script:
    """
    minimap2 -a \\ 
             -x asm20 \\ # optimized for mapping genome assemblies
             --score-N=0 \\ # make alignments neutral to N (no penalty)
             -N 0 \\ # no secondary alignments
             -p 0.8 \\ 
             -o ${ref_genome} ${selected_segments} > alignment.sam
    """

}
