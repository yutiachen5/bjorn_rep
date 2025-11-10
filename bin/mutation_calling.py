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

    return match["GFF_FEATURE"].to_list() if match.height > 0 else []

def get_codon_aa(ref_seq, alt_seq, pos): # pos: 1-based
    start = ((pos-1) // 3)*3
    ref_codon = ref_seq[start: start + 3]
    alt_codon = alt_seq[start: start + 3]

    ref_aa = Seq(ref_codon).translate(table=1)
    alt_aa = Seq(alt_codon).translate(table=1)

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

            # skip the seq with more than x% gaps or missing values
            # if (alt_seq.count("-")+alt_seq.count("N")) / length > 0.77:
            #     continue

            while i < length:
                if ref_seq[i] != alt_seq[i]:
                    gff_feature = get_gff_feature(gff, i+1)

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
                        
                        if gff_feature:
                            for g in gff_feature:
                                # ref_codon, alt_codon, ref_aa, alt_aa = get_codon_aa(ref_seq, alt_seq, i+1) 
                                outfile.write(f"{record.id}\t{args.region}\t{i}\t{ref_seq[i-1]}\t-{deleted}\t{g}\n") # use the preceding nuc as ref
                        else:
                            outfile.write(f"{record.id}\t{args.region}\t{i}\t{ref_seq[i-1]}\t-{deleted}\t""\n") 

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
                        
                        if gff_feature:
                            for g in gff_feature:
                                # ref_codon, alt_codon, ref_aa, alt_aa = get_codon_aa(ref_seq, alt_seq, i+1) 
                                outfile.write(f"{record.id}\t{args.region}\t{i+1}\t-\t+{inserted}\t{g}\n")
                        else:
                            outfile.write(f"{record.id}\t{args.region}\t{i+1}\t-\t+{inserted}\t""\n")

                        i = j + 1
                        continue

                    # SNP
                    else: 
                        if gff_feature:
                            for g in gff_feature:
                                # ref_codon, alt_codon, ref_aa, alt_aa = get_codon_aa(ref_seq, alt_seq, i+1) 
                                outfile.write(f"{record.id}\t{args.region}\t{i+1}\t{ref_seq[i]}\t{alt_seq[i]}\t{g}\n") 
                        else:
                            outfile.write(f"{record.id}\t{args.region}\t{i+1}\t{ref_seq[i]}\t{alt_seq[i]}\t""\n") 
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
#   -a /home/eleanor124/projects/bjorn_rep/output/PB2/alignment.fasta \
#   --gff /home/eleanor124/projects/bjorn_rep/data/PB2-DMS/PP755596.1.gff \
#   --region PB2 \
#   --ref_file /home/eleanor124/projects/bjorn_rep/data/Hu1-BA/ref_all.fasta \
#   --ref_id NC_045512.2

