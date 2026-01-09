#!/usr/bin/env python

import re
import argparse
import itertools
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

IUPAC = {
    "R": ["A","G"],
    "Y": ["C","T"],
    "S": ["G","C"],
    "W": ["A","T"],
    "K": ["G","T"],
    "M": ["A","C"],
    "B": ["C","G","T"],
    "D": ["A","G","T"],
    "H": ["A","C","T"],
    "V": ["A","C","G"],
    "N": ["A","C","G","T"],
    "-": ["-"],
}

def expand_base(b):
    if b in "ACGT":
        return [b]
    else:
        return IUPAC[b]
    
def translate(codon):
    try:
        return str(Seq("".join(codon)).translate(table=1))
    except:
        return "X"

def get_possible_aas(codon):
    possibilities = itertools.product(
        expand_base(codon[0]),
        expand_base(codon[1]),
        expand_base(codon[2]),
    )
    return {translate(p) for p in possibilities}

def consensus_aa(codon):
    aas = get_possible_aas(codon)
    return next(iter(aas)) if len(aas) == 1 else "X"

def extract_sid(ref_id, sra):
    # SC2
    if re.search(r"NC_045512\.2", ref_id, flags=re.IGNORECASE) is not None:
        if len(sra.split("/")) == 4:
            sid = sra.split("/")[2]
        elif sra.startswith("Consensus_"):
            sid = sra.split("_")[1]
        else:
            sid = sra
    # H5N1
    else:
        sid = sra.split("_")[1] # TODO: handle other viruses 
    return sid

def parse_ref_file(args):
    ref_records = list(SeqIO.parse(args.ref_file, "fasta"))
    matched = [rec for rec in ref_records if rec.id == args.ref_id]

    if len(matched) != 1:
        raise ValueError(f"Reference ID '{args.ref_id}' not found or multiple records found in {args.ref_file}")

    ref_seq = str(matched[0].seq).upper()
    
    return ref_seq, len(ref_records)

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
                gff.append({"local_start": start, "local_end": end, "GFF_FEATURE": gff_feature})

    gff_df = pd.DataFrame(gff)
    gff_df_grouped = (
        gff_df.groupby("GFF_FEATURE")
            .agg(global_start=("local_start", "min"),
                 global_end=("local_end", "max"))
            .reset_index()
    )

    gff = gff_df.merge(gff_df_grouped, on="GFF_FEATURE", how="left")

    return gff

def get_gff_feature(gff, pos):
    # GFF intervals are 1-based, inclusive
    # pos is 1-based index
    match = gff[(gff["local_start"] <= pos) & (pos <= gff["local_end"])]

    if len(match) > 0:
        gff_feature = match["GFF_FEATURE"].tolist()
        local_start = match["local_start"].tolist()
        local_end = match["local_end"].tolist()
        global_start = match["global_start"].tolist()
        global_end = match["global_end"].tolist()
    else:
        gff_feature, local_start, local_end, global_start, global_end = [], [], [], [], []

    return gff_feature, local_start, local_end, global_start, global_end

def get_codon_aa(ref_seq, alt_seq, pos, local_start, local_end, global_start, global_end): 
    # pos, start, end are all 1-based idx
    codon_start = (pos - local_start) // 3 * 3 + local_start - 1
    if codon_start + 2 > local_end:
        raise Exception ("codon index exceeds the CDS interval")

    ref_codon = ref_seq[codon_start: codon_start + 3]
    alt_codon = alt_seq[codon_start: codon_start + 3]

    ref_aa = consensus_aa(ref_codon)
    alt_aa = consensus_aa(alt_codon)

    if global_start != local_start:
        pos_aa = (pos - local_start) // 3 + 1 + (local_start - global_start) // 3 + 1 
    else:
        pos_aa = (pos - local_start) // 3 + 1

    return ref_codon, alt_codon, ref_aa, alt_aa, pos_aa

def mutation_calling(args):
    records = list(SeqIO.parse(args.a, "fasta"))
    ref_seq, n_ref = parse_ref_file(args)
    gff = parse_gff_file(args.gff)
    
    seq_len = [len(record.seq) for record in records]
    assert len(set(seq_len)) == 1, "Sequence lengths mismatch!"
    length = seq_len.pop()

    mut = []
    helper = []
    
    header = ["sra", "region", "pos", "ref", "alt", "GFF_FEATURE", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "pos_aa"]
    helper_header = ["sra", "pos", "alt"]

    for record in records[n_ref: ]: # refs are put at the top of this file
        i = 0
        alt_seq = str(record.seq).upper()
        sid = extract_sid(args.ref_id, record.id)

        while i < length:
            if ref_seq[i] != alt_seq[i]:
                gff_feature, local_start, local_end, global_start, global_end = get_gff_feature(gff, i+1) 

                # DEL
                if alt_seq[i] == "-": 
                    j = i
                    while j + 1 < length and alt_seq[j+1] == "-":
                        j += 1
                    deleted = ref_seq[i:j+1]

                    if i == 0 or j == length - 1: # leading or trailing del
                        helper.append([sid, i, "-"+deleted])
                        i = j + 1
                        continue
                    
                    mut.append([sid, args.region, i, ref_seq[i-1], "-"+deleted, "", "", "", "", "", np.nan])
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
                    mut.append([sid, args.region, i, ref_seq[i-1], "+"+inserted, "", "", "", "", "", np.nan])
                    i = j + 1
                    continue

                # SNP
                else: 
                    if gff_feature:
                        flag = False
                        for k in range(len(gff_feature)):
                            ref_codon, alt_codon, ref_aa, alt_aa, pos_aa = get_codon_aa(ref_seq, alt_seq, i+1, local_start[k], local_end[k], global_start[k], global_end[k]) 
                            # if ref_codon contains N or gap, skip 
                            if "N" in ref_codon or "-" in ref_codon:
                                continue
                            # synonymous mutation and unknown aa
                            if ref_aa == alt_aa: 
                                if flag == False:
                                    mut.append([sid, args.region, i+1, ref_seq[i], alt_seq[i], "", "", "", "", "", np.nan])
                                    flag = True
                                else:
                                    continue
                            else:
                                mut.append([sid, args.region, i+1, ref_seq[i], alt_seq[i], gff_feature[k]+'_'+args.region, ref_codon, alt_codon, ref_aa, alt_aa, pos_aa])
                    else:
                        mut.append([sid, args.region, i+1, ref_seq[i], alt_seq[i], "", "", "", "", "", np.nan])
            i += 1

    helper = pd.DataFrame(helper, columns=helper_header) 
    helper = helper.astype({
        "sra": "string",
        "pos": "Int64",
        "alt": "string",
    })      
    helper.to_csv("del_helper.tsv", sep="\t", index=False, header=helper_header)

    mut = pd.DataFrame(mut, columns=header)
    mut = mut.astype({
        "sra": "string",
        "region": "string",
        "pos": "Int64",
        "ref": "string",
        "alt": "string",
        "GFF_FEATURE": "string",
        "ref_codon": "string",
        "alt_codon": "string",
        "ref_aa": "string",
        "alt_aa": "string",
        "pos_aa": "Int64",
    })
    mut.to_csv(args.o, sep="\t", index=False, header=header)


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

# python /home/yutianc/bjorn_rep/bin/mutation_calling.py \
#   -a /home/yutianc/bjorn_rep/output/Hu1/mm/alignment.fasta \
#   --gff /home/yutianc/bjorn_rep/data/Hu1-BA/NC_045512.2.gff \
#   --region NC_045512.2 \
#   --ref_file /home/yutianc/bjorn_rep/data/Hu1-BA/alignment.fasta \
#   --ref_id NC_045512.2
