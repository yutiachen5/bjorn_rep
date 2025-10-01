#!/usr/bin/env python

import os
import pandas as pd
import argparse
import polars as pl
from Bio import SeqIO

def extract_variants(args):
    # background genome should be put at the first line of this file
    records = list(SeqIO.parse(args["a"], "fasta"))

    ids = [args["ref"], args["query"]]
    rec_dic = {record.id: record.seq for record in records if record.id in ids}
    seq_len = set([len(record.seq) for record in records if record.id in ids])

    assert len(seq_len) == 1, "Alignment length mismatch" # currently only for the case when all genome havce the same length

    n = seq_len.pop()
    variants = [] # {"pos":int, "gid1":str, ...}

    for i in range(n):
        nuc = {sid: seq[i].upper() for sid, seq in rec_dic.items()} 
        if len(set(nuc.values())) > 1: # mismatch
            row = {"pos": i+1} # 1-based idx
            row.update(nuc)
            variants.append(row)

    return pl.DataFrame(variants)

def derive_mutations(variants, args):
    mutation = pl.read_csv(args["m"], separator="\t", has_header=True)
    dftmp = variants.select(
        ["pos", args["ref"], args["query"]]
    ).filter(
        pl.col(args["ref"]) != pl.col(args["query"])
    )

    merged = dftmp.join(mutation, on='pos', how='right')

    # if bg_genome.alt = alt_genome.ref: no mutation, rm the line
    # else: different mutation, change ref nuc
    merged = (merged.with_columns(
        pl.when(pl.col("alt") == pl.col(args["query"]))
            .then(None)
            .otherwise(pl.col("alt"))
            .alias("alt"),

        pl.when(pl.col("alt") != pl.col(args["query"]))
            .then(pl.col(args["query"]))
            .otherwise(pl.col("ref"))
            .alias("ref")
    )
    .filter(pl.col("alt").is_not_null())
    .drop([args["ref"], args["query"]])
    .sort("pos")
    )

    merged.write_csv(args["query"]+"_"+args["o"], separator="\t")


def main():
    parser = argparse.ArgumentParser(description="Extract differences between Hu-1 and other references to infer mutations on other refences.")
    parser.add_argument("-m", help="Mutation file path in TSV format.", required=True)
    parser.add_argument("-a", help="Alignment file from mafft.", required=True)
    parser.add_argument("-o", help="Output diractory for new mutation file.")
    parser.add_argument("--ref", help="ID of refence genome", required=True)
    parser.add_argument("--query", help="ID of query genom", required=True)
    args = vars(parser.parse_args())

    variants = extract_variants(args)
    derive_mutations(variants, args)

    
if __name__ == "__main__":
    main()