#!/usr/bin/env python

import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import polars as pl

def parse_gff_file(gff_path):
    gff = []
    with open(gff_path) as gfile:
        for line in gfile:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                continue
            if cols[2] == 'CDS': # cds only
                start, end = int(cols[3]), int(cols[4])
                gff_feature = re.search(r"(?<=cds-)[^;\s]+", cols[-1]).group()
                gff.append({"start": start, "end": end, "GFF_FEATURE": gff_feature})

    return pl.DataFrame(gff)

def get_gff_feature(gff, pos): # pos: 1-based
    match = gff.filter((pl.col("start") <= pos) & (pos <= pl.col("end"))) # GFF intervals are 1-based, inclusive

    if match.height > 0:
        gff_feature = match["GFF_FEATURE"].to_list()
        start = match["start"].to_list()
        end = match["end"].to_list()
    else:
        gff_feature, start, end = [], [], []

    return gff_feature, start, end

def get_codon_aa(ref_seq, alt_seq, pos, start, end): # pos: 1-based
    codon_start = (pos-start)//3*3 + start - 1
    if codon_start + 3 > end:
        raise Exception ("codon index exceeds the CDS interval")

    ref_codon = ref_seq[codon_start: codon_start + 3]
    alt_codon = alt_seq[codon_start: codon_start + 3]

    ref_aa = Seq(ref_codon).translate(table=1)
    try:
        alt_aa = Seq(alt_codon).translate(table=1)
    except:
        alt_aa = ""

    return ref_codon, alt_codon, ref_aa, alt_aa

def parse_ref_file(args):
    ref_records = list(SeqIO.parse(args.ref_file, "fasta"))
    matched = [rec for rec in ref_records if rec.id == args.ref_id]

    if len(matched) != 1:
        raise ValueError(f"Reference ID '{args.ref_id}' not found or multiple records found in {args.ref_file}")

    ref_seq = str(matched[0].seq)
    
    return ref_seq, len(ref_records)


def mutation_calling(args):
    records = list(SeqIO.parse(args.a, "fasta"))
    ref_seq, n_ref = parse_ref_file(args)

    gff = parse_gff_file(args.gff)
    
    seq_len = [len(record.seq) for record in records]
    assert len(set(seq_len)) == 1, "Sequence lengths mismatch!"

    length = seq_len.pop()
    

    with open(args.o, "w", newline="") as outfile, open("del_helper.tsv", "w", newline="") as helper:
        header = ["sra", "region", "pos", "ref", "alt", "GFF_FEATURE", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "pos_aa"]
        outfile.write("\t".join(header) + "\n")
        helper.write("\t".join(header) + "\n")

        for record in records[n_ref: ]: # refs are put at the top of this file
            i = 0
            alt_seq = str(record.seq).upper()

            while i < length:
                if ref_seq[i] != alt_seq[i]:
                    gff_feature, start, end = get_gff_feature(gff, i+1) # only do gff feature for SNP???

                    # DEL
                    if alt_seq[i] == "-": 
                        j = i
                        while j + 1 < length and alt_seq[j+1] == "-":
                            j += 1
                        deleted = ref_seq[i:j+1]

                        if i == 0 or j == length - 1: # leading or trailing del
                            helper.write(f"{record.id}\t{args.region}\t{i}\t{'NA'}\t-{deleted}\t""\n")
                            i = j + 1
                            continue
                        
                        outfile.write(f"{record.id}\t{args.region}\t{i}\t{ref_seq[i-1]}\t-{deleted}\t""\t""\t""\t""\t""\n") 

                        i = j + 1
                        continue

                    # INS
                    elif ref_seq[i] == "-": 
                        j = i
                        while j + 1 < length and ref_seq[j+1] == "-":
                            j += 1
                        inserted = alt_seq[i:j+1]

                        # Skip leading or trailing insertions -- keep???
                        # if i == 0 or j == length - 1:
                        #     i = j + 1
                        #     continue

                        outfile.write(f"{record.id}\t{args.region}\t{i+1}\t-\t+{inserted}\t""\t""\t""\t""\t""\n")

                        i = j + 1
                        continue

                    # SNP
                    else: 
                        if gff_feature:
                            for k in range(len(gff_feature)):
                                ref_codon, alt_codon, ref_aa, alt_aa = get_codon_aa(ref_seq, alt_seq, i+1, start[k], end[k]) 
                                outfile.write(f"{record.id}\t{args.region}\t{i+1}\t{ref_seq[i]}\t{alt_seq[i]}\t{gff_feature[k]}\t{ref_codon}\t{alt_codon}\t{ref_aa}\t{alt_aa}\n") 
                        else:
                            outfile.write(f"{record.id}\t{args.region}\t{i+1}\t{ref_seq[i]}\t{alt_seq[i]}\t""\t""\t""\t""\t""\n") 
                i += 1
    


def main():
    parser = argparse.ArgumentParser(description="Extract differences between Hu-1 and other references to infer mutations on other refences.")
    parser.add_argument("-a", help="Alignment file path from gofasta in FASTA format, including bakground genome and query genomes.", required=True)
    parser.add_argument("--gff", help="GFF file to extract gene features and sra.", required=True)
    parser.add_argument("-o", help="Output name for new mutation file in TSV format.", default="mutations.tsv")
    parser.add_argument("--region", help="Region of reference genome.", required=True, type=str)
    parser.add_argument("--ref_file", help="FASTA file including background and query genome in reference.", required=True, type=str)
    parser.add_argument("--ref_id", help="ID of background reference genome.", required=True, type=str)


    args = parser.parse_args()

    mutation_calling(args)
    

    
if __name__ == '__main__':
    main()

# python /home/eleanor124/projects/bjorn_rep/bin/mutation_calling.py \
#   -a /home/eleanor124/projects/bjorn_rep/testing/5YEQUKFCG/alignment.fasta \
#   --gff /home/eleanor124/projects/bjorn_rep/data/Hu1-BA/NC_045512.2.gff \
#   --region Hu1 \
#   --ref_file /home/eleanor124/projects/bjorn_rep/data/Hu1-BA/NC_045512.2_BA.1.fasta \
#   --ref_id NC_045512.2

