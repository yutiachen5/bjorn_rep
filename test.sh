#!/usr/bin/env bash

REF="$1"
QUERY="$2"
OUT_PREFIX="$3"

minimap2 -a -x map-ont --score-N=0 --secondary=no "$REF" "$QUERY" > "${OUT_PREFIX}.sam" 

gofasta sam toMultiAlign \
    -s "${OUT_PREFIX}.sam" \
    -r "${REF}" \
    -o "${OUT_PREFIX}.fasta"