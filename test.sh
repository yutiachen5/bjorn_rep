#!/usr/bin/env bash

REF="$1"
QUERY="$2"
GFF="$3"
OUT_PREFIX="$4"

mkdir -p ./testing/${OUT_PREFIX}

cd ./testing/${OUT_PREFIX}

minimap2 -a -x asm20 --score-N=0 --sam-hit-only --secondary=no "$REF" "$QUERY" > "alignment.sam" 

gofasta sam toMultiAlign \
    -s "alignment.sam" \
    -r "${REF}" \
    -o "alignment.fasta"

gofasta variants --msa "alignment.fasta" \
                    --reference CY018884.1_cds_ABM21959.1_1 \
                    -a "${GFF}" \
                    --append-snps \
                    -o "aa_changes.csv"
