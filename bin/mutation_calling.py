#!/usr/bin/env python

import re
import argparse
from Bio import SeqIO
import pandas as pd
import polars as pl

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
                start, end = int(cols[3]), int(cols[4])
                gff_feature = re.search(r"YP_\d+\.\d+", cols[-1]).group() 
                gff.append({"start": start, "end": end, "GFF_FEATURE": gff_feature})

    return pl.DataFrame(gff)


def mutation_calling(args):
    records = list(SeqIO.parse(args.a, "fasta"))
    
    seq_len = [len(record.seq) for record in records]
    assert len(set(seq_len)) == 1, "Sequence lengths mismatch!"

    length = seq_len.pop()

    mutations = []
    for record in records[1:]:
        # mut_type: "<type>:<start>:<length>"

        i = 0

        ref_seq = str(records[0].seq)
        alt_seq = str(record.seq)

        while i < length:
            if ref_seq[i] != alt_seq[i]:
                if alt_seq[i] == "-": # del
                    j = i
                    while j + 1 < length and alt_seq[j+1] == "-":
                        j += 1
                    deleted = ref_seq[i:j+1]
                    mutations.append({"sra": record.id, "region": args.region, "pos": i+1, "ref": ref_seq[i], "alt": "-" + deleted, "mut_len": len(deleted)})

                    i = j+1
                    continue
                elif ref_seq[i] == "-": # ins
                    j = i
                    while j + 1 < length and ref_seq[j+1] == "-":
                        j += 1
                    inserted = alt_seq[i:j+1]
                    mutations.append({"sra": record.id, "region": args.region, "pos": i+1, "ref": ref_seq[i], "alt": "+" + inserted, "mut_len": len(inserted)})

                    i = j+1
                    continue
                else: # snp
                    mutations.append({"sra": record.id, "region": args.region, "pos": i+1, "ref": ref_seq[i], "alt": alt_seq[i], "mut_len": 1})
            i += 1

    schema = {
        "sra": pl.Utf8,
        "region": pl.Utf8,
        "pos": pl.Int64,
        "ref": pl.Utf8,
        "alt": pl.Utf8,
        "mut_len": pl.Int64,
    }
    mutations = pl.DataFrame(mutations, schema=schema)
    gff = parse_gff(args.gff)

    mutations = (
        mutations.sort("pos")
            .join_asof(
                gff.sort("start"),
                left_on="pos",
                right_on="start",
                strategy="backward" # pos in variants_sample was matched with the closet start position from gff, where pos>start
            )
            .with_columns(
                pl.when(pl.col("pos") < pl.col("end"))
                .then(pl.col("GFF_FEATURE"))
                .otherwise(None)
                .alias("GFF_FEATURE")
            )
            .with_columns([
                pl.lit(None).alias("ref_codon"),
                pl.lit(None).alias("alt_codon"),
                pl.lit(None).alias("ref_aa"),
                pl.lit(None).alias("alt_aa"),
                pl.lit(None).alias("pos_aa")
            ])
            .drop(["start", "end"])
    )
    
    mutations.write_csv(args.o, separator="\t", include_header=True)


def main():
    parser = argparse.ArgumentParser(description="Extract differences between Hu-1 and other references to infer mutations on other refences.")
    parser.add_argument("-a", help="Alignment file path from gofasta in FASTA format, including bakground genome and query genomes.", required=True)
    parser.add_argument("--gff", help="GFF file to extract gene features and sra.", required=True)
    parser.add_argument("-o", help="Output name for new mutation file in TSV format.", default="mutations.tsv")
    parser.add_argument("--region", help="Region of reference genome.", required=True, type=str)

    args = parser.parse_args()

    mutation_calling(args)

    
if __name__ == '__main__':
    main()

# python /home/eleanor124/projects/bjorn_rep/bin/mutation_calling.py \
#   -a /home/eleanor124/projects/bjorn_rep/data/testing/aligned.fasta \
#   --gff /home/eleanor124/projects/bjorn_rep/data/NC_045512.2.gff \
#   --region NC_045512.2 \
#   -o /home/eleanor124/projects/bjorn_rep/data/testing/mutations.tsv