#!/usr/bin/env bash

REF="$1"
QUERY="$2"
OUT_PREFIX="$3"

minimap2 -a -t 1 -x asm20 --score-N=0 --secondary=no "$REF" "$QUERY" > "${OUT_PREFIX}.sam" 

gofasta sam toMultiAlign \
    -s "${OUT_PREFIX}.sam" \
    -r "${REF}" \
    -o "${OUT_PREFIX}.fasta"

# gofasta_converter.py --csv_path ${aa_changes_csv} \\
#                     --reference ${params.ref_file} \\
#                     --region ${params.region} \\                    
#                     --output mutations.tsv