#!/usr/bin/env python

import pandas as pd
import argparse
import polars as pl
from Bio import SeqIO

pl.Config.set_tbl_cols(-1)      

def extract_variants(args):
    # background genome should be put at the first line of this file
    records = list(SeqIO.parse(args["a"], "fasta"))

    ids = [args["bg"], args["query"]]
    rec_dic = {record.id: record.seq for record in records if record.id in ids}
    seq_len = set([len(record.seq) for record in records if record.id in ids])

    assert len(seq_len) == 1, "Alignment length mismatch" # currently only for the case when all genome havce the same length

    n = seq_len.pop()
    variants = [] # {"pos":int, "gid1":str, ...}

    for i in range(n):
        nuc = {gid: seq[i].upper() for gid, seq in rec_dic.items()} # {"NC_...":"ATCG..", "NC_...BA1":"ATTG"}
        if nuc[args["bg"]] != nuc[args["query"]]:
            row = {"pos": i+1} # 1-based idx
            row.update(nuc)
            variants.append(row)

    return pl.DataFrame(variants)

def derive_mutations(variants, args):
    mutation = pl.read_csv(args["m"], separator="\t", has_header=True)
    samples = mutation.select(["sra", "region"].unique()) # unique samples

    B = args["bg"] # bg.ref
    Q = args["query"] # q.ref

    both = (pl.col("sra").is_not_null()) & (pl.col(B).is_not_null())
    left_only = (pl.col("sra").is_not_null()) & (pl.col(B).is_null())
    right_only = (pl.col("sra").is_null()) & (pl.col(B).is_not_null())

    merged = mutation.join(variants, on="pos", how="full") # left: mutation, right: variants

    assert merged.filter(both).select((pl.col(B) == pl.col("ref")).all()).item(), "mismatch between ref nucs!"
    assert merged.select((pl.col("ref") != pl.col("alt")).all()).item(), "ref shouldn't match alt in mutation file!"

    merged = merged.rename({"ref": "old_ref", "alt": "old_alt"})   
                
    # bg.ref != q.ref already guaranteed
    # only consider the record in the mutation file
    # case1: both & bg.alt == q.ref, rm
    # case2: both & bg.alt != q.ref & bg.alt != bg.ref, update: ref = q.ref
    # case3: variants only (bg.alt == bg.ref != r.ref - mutation DNE), add: ref = q.ref, alt = bg.ref 
    # case4: mutation only (bg.ref == q.ref), keep


    res = (
        merged.with_columns(
            # case1
            pl.when(both & (pl.col("old_alt") == pl.col(Q)))
                .then(None)

            # case2
            .when(both & (pl.col("old_alt") != pl.col(Q)))
                .then(pl.struct([pl.col(Q).alias("ref"),
                                 pl.col("old_alt").alias("alt")]))
            # case3
            .when(right_only)
                .then(pl.struct([pl.col(Q).alias("ref"),
                                 pl.col(B).alias("alt")]))
            # case4
            .when(left_only)
                .then(pl.struct([pl.col("old_ref").alias("ref"),
                                 pl.col("old_alt").alias("alt")]))
            .otherwise(None)  
            .alias("case")
        )
        .filter(pl.col("case").is_not_null())
        .coalesce([pl.col("pos"), pl.col("pos_right")]).alias("pos")
        .drop([B, Q, "old_ref", "old_alt", "pos_right"])      
        .unnest("case")          # expand struct into ref/alt columns
        .sort("pos")
    )

    assert merged.filter(both).select((pl.col("pos") == pl.col("pos_right")).all()).item(), "mismatch between pos!"

    sample_appended = res.filter(right_only).join()
    res = res.join(samples, how="corss")


    print(res.head(5))
    print(res.height)
    # merged.write_csv(Q+"_"+args["o"], separator="\t")

    # res.write_csv("tmp.tsv", separator="\t")


def main():
    parser = argparse.ArgumentParser(description="Extract differences between Hu-1 and other references to infer mutations on other refences.")
    parser.add_argument("-m", help="Mutation file path in TSV format.", required=True)
    parser.add_argument("-a", help="Alignment file path from mafft.", required=True)
    parser.add_argument("-o", help="Output diractory for new mutation file.")
    parser.add_argument("--bg", help="ID of background genome", required=True)
    parser.add_argument("--query", help="ID of query genom", required=True)
    args = vars(parser.parse_args())

    variants = extract_variants(args)
    derive_mutations(variants, args)

    
if __name__ == "__main__":
    main()

# for testing locally:
# python /home/eleanor124/projects/bjorn_rep/bin/translate_mutations.py -m /home/eleanor124/projects/bjorn_rep/output/mutations.tsv -a /home/eleanor124/projects/bjorn_rep/mafft_oneline.fasta --query 'NC_045512.2_BA.1' --bg 'NC_045512.2'