#!/usr/bin/env python

import re
import argparse
import polars as pl
from Bio import SeqIO

pl.Config.set_tbl_cols(-1)   

def extract_ref_variants(args): 
    records = list(SeqIO.parse(args.a, "fasta"))

    seq_len = set([len(record.seq) for record in records])
    assert len(seq_len) == 1, "Sequence length mismatch!"

    ids = [args.bg, args.query]
    rec_dic = {record.id: record.seq for record in records if record.id in ids}

    references = []
    for i in range(seq_len.pop()): 
        if rec_dic[args.bg][i] != rec_dic[args.query][i]:
            references.append({"pos": i+1, "bg.ref":rec_dic[args.bg][i], "q.ref":rec_dic[args.query][i]}) # 1-based index

    return pl.DataFrame(references)

def parse_gff(gff_path):
    gff = []
    with open(gff_path) as gfile:
        for line in gfile:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                continue
            if cols[2] == 'CDS': # cds only, ID=id-YP_009724389.1:5325..5925
                region = cols[0]
                start, end = int(cols[3]), int(cols[4])
                gff_feature = re.search(r"YP_\d+\.\d+", cols[-1]).group() 
                gff.append({"region":region, "start": start, "end": end, "GFF_FEATURE": gff_feature})

    return pl.DataFrame(gff)

def translate_mutations(args):
    references = extract_ref_variants(args)
    mutation = pl.read_csv(args.m, separator="\t", has_header=True)
    sample_id = mutation.select(pl.col("sra").unique())
    gff = parse_gff(args.gff)
    

    both = (pl.col("sra").is_not_null()) & (pl.col("bg.ref").is_not_null())
    left_only = (pl.col("sra").is_not_null()) & (pl.col("bg.ref").is_null())
    right_only = (pl.col("sra").is_null()) & (pl.col("bg.ref").is_not_null())

    merged = mutation.join(references, on="pos", how="full") # left: mutation, right: references

    print(merged.filter(both & (pl.col("bg.ref") != pl.col("ref"))))
    assert (
        merged.filter(both & pl.col("alt").str.len_chars() == 1)
            .select((pl.col("bg.ref") == pl.col("ref")).all())
            .item()
    ), "Mismatch between ref nucs for SNPs!"
    assert (merged.select((pl.col("ref") != pl.col("alt")).all()).item()), "ref shouldn't match alt in mutation file!"

    merged = merged.rename({"ref": "old_ref", "alt": "old_alt"})   
                
    # case1: both & old_alt == q.ref, rm
    # case2: both & old_alt != q.ref
        # case2a - snp: update: ref = q.ref, alt = old_alt
        # case2b - ins: update: ref = q.ref, alt = old_alt
        # case2c - del: ref = q.ref, alt = -q.ref
    # case3: mut on new ref only (old_alt == bg.ref != q.ref), add: ref = q.ref, alt = bg.ref 
    # case4: mut on old ref only (bg.ref == q.ref), keep: ref = old_ref, alt = old_alt

    merged = (
        merged.with_columns(
            # case1
            pl.when(both & (pl.col("old_alt") == pl.col("q.ref")))
                .then(None)

            # case2a and 2b - ins and snp
            .when(both & (pl.col("old_alt") != pl.col("q.ref")) & (pl.col("old_alt").str.starts_with("+") | pl.col("old_alt").str.len_chars() == 1))
                .then(pl.struct([pl.col("q.ref").alias("ref"),
                                 pl.col("old_alt").alias("alt")]))
            # case2c - del
            .when(both & (pl.col("old_alt").str.starts_with("-")))
                .then(pl.struct([pl.col("q.ref").alias("ref"),
                                 ("-" + pl.col("q.ref")).alias("alt")]))
            # case3
            .when(right_only)
                .then(pl.struct([pl.col("q.ref").alias("ref"),
                                 pl.col("bg.ref").alias("alt")]))
            # case4
            .when(left_only)
                .then(pl.struct([pl.col("old_ref").alias("ref"),
                                 pl.col("old_alt").alias("alt")]))
            .otherwise(None)  
            .alias("case"),
            pl.coalesce([pl.col("pos"), pl.col("pos_right")]).alias("pos")
        )
        .filter(pl.col("case").is_not_null())
        .drop(["q.ref", "old_ref", "old_alt", "pos_right"])      
        .unnest("case")          
        .sort("pos")
    )

    # remove the sample if its mut pos is inside the region of del
    deletion = (
        merged.filter(pl.col("alt").str.starts_with("-"))
           .with_columns((pl.col("pos") + pl.col("mut_len") - 1).alias("end"))
           .select(["sra", "pos", "end"])
           .sort(["sra", "pos"])
    )

    right = merged.filter(right_only)
    left_both = merged.filter(~right_only)

    # get the sample info for the variants coming from new ref genome
    right = right.join(sample_id, how="cross")
    right = (
        right.sort("pos")
            .join_asof(
                gff.sort("start"),
                left_on="pos",
                right_on="start",
                strategy="backward" # pos in variants_sample was matched with the closet start position from gff, where pos>start
            )
            .with_columns(
                pl.when(pl.col("pos") <= pl.col("end"))
                    .then(pl.col("GFF_FEATURE_right"))
                    .otherwise(pl.col("GFF_FEATURE"))
                    .alias("GFF_FEATURE")
            )
            .with_columns([
                pl.coalesce([pl.col("sra"), pl.col("sra_right")]).alias("sra"),
                pl.coalesce([pl.col("region"), pl.col("region_right")]).alias("region")
            ])
            .drop(["sra_right", "start", "end", "GFF_FEATURE_right", "region_right"])
            .sort(["sra", "pos"])
    )

    right = (
        right.join_asof(
            deletion.sort(["sra", "pos"]),
            left_on="pos",
            right_on="pos",
            by="sra", # join per sample
            strategy="backward",
        )
        .filter(~(pl.col("pos") <= pl.col("end")))  # keep only those NOT in deletion
        .drop(["end"])
    )

    final = (
        pl.concat([right, left_both])
                .drop(["bg.ref"]) # todo: drop mut_len in output files??
                .sort("pos")
    )

    # final.write_csv("tmp.tsv", separator="\t", include_header=True)
    final.write_csv(args.query+"_"+args.o, separator="\t", include_header=True)


def main():
    parser = argparse.ArgumentParser(description="Extract differences between Hu-1 and other references to infer mutations on other refences.")
    parser.add_argument("-m", help="Mutation file path in TSV format.", required=True)
    parser.add_argument("-a", help="Alignment file path from gofasta in FASTA format, including bakground genome and query genomes.", required=True)
    parser.add_argument("--gff", help="GFF file to extract gene features and sra.", required=True)
    parser.add_argument("-o", help="Output name for new mutation file in TSV format.", default="mutations.tsv")
    parser.add_argument("--bg", help="ID of background genome", required=True, type=str)
    parser.add_argument("--query", help="ID of query genom", required=True, type=str)
    args = parser.parse_args()

    translate_mutations(args)

    
if __name__ == '__main__':
    main()

# cmd to test locally:
# python /home/eleanor124/projects/bjorn_rep/bin/translate_mutations.py \
#   -m /home/eleanor124/projects/bjorn_rep/data/testing/mutations.tsv \
#   -a /home/eleanor124/projects/bjorn_rep/data/testing/HU1_BA1_BA2_final.fa \
#   --gff /home/eleanor124/projects/bjorn_rep/data/NC_045512.2.gff \
#   --bg NC_045512.2 \
#   --query NC_045512.2_BA.1 