#!/usr/bin/env python

import re
import argparse
import itertools
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

pd.set_option("display.max_columns", None)


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
    codon = codon.upper()
    possibilities = itertools.product(
        expand_base(codon[0]),
        expand_base(codon[1]),
        expand_base(codon[2]),
    )
    return {translate(p) for p in possibilities}

def consensus_aa(codon):
    aas = get_possible_aas(codon)
    return next(iter(aas)) if len(aas) == 1 else "X"

def parse_alignment_file(args): 
    records = list(SeqIO.parse(args.a, "fasta"))

    seq_len = set([len(record.seq) for record in records])
    assert len(seq_len) == 1, "Sequence length mismatch!"
    length = seq_len.pop()

    ref_records = records[:args.n_ref] # refs are put at the top of alignment file
    sample_records = records[args.n_ref: ]

    ref_dic = {record.id: str(record.seq) for record in ref_records if record.id in [args.bg, args.query]}
    sample_dic = {record.id: str(record.seq) for record in sample_records}

    references = []
    for i in range(length): 
        if ref_dic[args.bg][i] != ref_dic[args.query][i]:
            references.append({"pos": i+1, "bg.ref":ref_dic[args.bg][i], "q.ref":ref_dic[args.query][i]}) # 1-based index

    return pd.DataFrame(references), ref_dic, sample_dic

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

def extract_del(mutation, helper):
    # region covered by normal del
    deletion = mutation[(mutation["alt"].str.startswith("-"))].copy()
    deletion["start"] = deletion["pos"] + 1
    deletion["end"] = deletion["pos"] + deletion["alt"].apply(lambda x: len(x)-1)

    # region covered by leading or trailing del, start and end are 1-based inclusive
    helper["start"] = helper["pos"] + 1
    helper["end"] = helper["pos"] + helper["alt"].apply(lambda x: len(x)-1)

    deletion = pd.concat([deletion, helper], axis=0)

    return deletion

def match_del(pos, s, deletion):
    match = deletion[(deletion["start"] <= pos) & (pos <= deletion["end"]) & (deletion["sra"] == s)]

    return True if not match.empty else False

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

    ref_codon = ref_seq[codon_start: codon_start + 3].upper()
    alt_codon = alt_seq[codon_start: codon_start + 3].upper()

    ref_aa = consensus_aa(ref_codon)
    alt_aa = consensus_aa(alt_codon)

    if global_start != local_start:
        pos_aa = (pos - local_start) // 3 + 1 + (local_start - global_start) // 3 + 1 
    else:
        pos_aa = (pos - local_start) // 3 + 1

    return ref_codon, alt_codon, ref_aa, alt_aa, pos_aa

def translate_mutations(args):
    mutation = pd.read_csv(args.m, sep="\t", header=0)
    helper = pd.read_csv(args.d, sep="\t", header=0)
    deletion = extract_del(mutation, helper)

    ref_variants, ref_dic, sample_dic = parse_alignment_file(args) # ref_variants: {"pos": 2, "Hu1": "A", "BA1": "T"}
    gff = parse_gff_file(args.gff)
    bg_seq, q_seq = ref_dic[args.bg], ref_dic[args.query]

    header = ["sra", "region", "pos", "ref", "alt", "GFF_FEATURE", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "pos_aa"]

    # case1: bg.ref == q.ref
        # case1a: mut_type == SNP, update ref_codon and ref_aa
        # case1b: mut_type == INS, appended all
        # case1c: mut_type == DEL, update ref seq
    # case2: bg.ref != q.ref
        # case2a: bg.ref != alt == q.ref, skip
        # case2b: bg.ref != alt != q.ref (alt can be del), update ref = q.ref, alt = alt
        # case2c: bg.ref == alt != q.ref, add ref = q.ref, alt = bg.ref

    pos_diff = ref_variants["pos"].to_list()
    pos_same = mutation[~mutation["pos"].isin(pos_diff)]

    # case1a
    pos_same_snp = pos_same[~(pos_same["alt"].str.startswith("-") | pos_same["alt"].str.startswith("+"))].copy()
    pos_same_snp["_tmp"] = 1
    gff["_tmp"] = 1

    pos_same_snp = pos_same_snp.merge(gff, on="_tmp", how="left", suffixes=("_drop", "")).drop(columns="_tmp")
    pos_same_snp = pos_same_snp[
        (pos_same_snp["local_start"] <= pos_same_snp["pos"]) &
        (pos_same_snp["pos"] <= pos_same_snp["local_end"])
    ]

    pos_same_snp = pos_same_snp[
        pos_same_snp["local_start"].isna()  |  (
            (pos_same_snp["pos"] >= pos_same_snp["local_start"]) &
            (pos_same_snp["pos"] <= pos_same_snp["local_end"])
        )
    ]#.reset_index(drop=True)

    pos_same_snp["codon_start"] = np.nan
    pos_same_snp["ref_codon"] = np.nan
    pos_same_snp["ref_aa"] = np.nan

    mask = pos_same_snp["local_start"].notna()
    pos_same_snp.loc[mask, "codon_start"] = (
        (pos_same_snp.loc[mask, "pos"] - pos_same_snp.loc[mask, "local_start"]) // 3 * 3+ pos_same_snp.loc[mask, "local_start"] - 1
    )
    pos_same_snp.loc[mask, "ref_codon"] = pos_same_snp.loc[mask].apply(
        lambda row: q_seq[int(row["codon_start"]): int(row["codon_start"]) + 3],
        axis=1
    )
    pos_same_snp.loc[mask, "alt_codon"] = pos_same_snp.loc[mask].apply(
        lambda row: sample_dic[row["sra"]][int(row["codon_start"]): int(row["codon_start"]) + 3],
        axis=1
    )

    # print(pos_same_snp[(pos_same_snp["sra"] == "Consensus_SRR24839088_PB2_cns_threshold_0.75_quality_20") & (pos_same_snp["pos"] == 1947)])
    pos_same_snp.loc[mask, "ref_aa"] = pos_same_snp.loc[mask, "ref_codon"].apply(consensus_aa)
    pos_same_snp.loc[mask, "alt_aa"] = pos_same_snp.loc[mask, "alt_codon"].apply(consensus_aa)

    pos_same_snp.loc[mask, "pos_aa"] = pos_same_snp.loc[mask].apply(
        lambda row: ((row["pos"] - row["local_start"]) // 3 + 1)
                    if row["global_start"] == row["local_start"]
                    else ((row["pos"] - row["local_start"]) // 3 + 1 +
                        (row["local_start"] - row["global_start"]) // 3 + 1),
        axis=1
    )

    # print(pos_same_snp[(pos_same_snp["sra"] == "Consensus_SRR24839088_PB2_cns_threshold_0.75_quality_20") & (pos_same_snp["pos"] == 1947)])

    pos_same_snp.loc[pos_same_snp["ref_aa"] == pos_same_snp["alt_aa"], ["GFF_FEATURE", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "pos_aa"]] = np.nan
    pos_same_snp = pos_same_snp.drop(columns=["codon_start", "local_start", "local_end", "global_start", "global_end", "GFF_FEATURE_drop"])

    pos_same_snp["pos_aa"] = pos_same_snp["pos_aa"].astype("Int64")
    pos_same_snp["pos_aa"] = pos_same_snp["pos_aa"].astype("Int64")

    # case1b
    pos_same_ins = pos_same[pos_same["alt"].str.startswith("+")]

    # case1c
    pos_same_del = pos_same[pos_same["alt"].str.startswith("-")].copy().reset_index(drop=True)
    if not pos_same_del.empty:
        pos_same_del["alt"] = pos_same_del.apply(lambda row: "-" + q_seq[row["pos"]: row["pos"] + len(row["alt"]) - 1], axis=1)

    pos_same_mut = pd.concat([pos_same_ins, pos_same_del, pos_same_snp], axis=0, ignore_index=True)


    pos_diff_mut = []
    for i in range(len(pos_diff)):
        pos = pos_diff[i] # 1-based index
        gff_feature, local_start, local_end, global_start, global_end = get_gff_feature(gff, pos)
        tmp = mutation[mutation["pos"] == pos]

        for s in list(sample_dic.keys()):
            del_match = match_del(pos, s, deletion)
            if del_match:
                continue
            if s in list(tmp["sra"]): # alt != bg.ref
                alt = tmp.loc[tmp["sra"] == s, "alt"].iloc[0]
                # case2a
                if alt == q_seq[pos-1]:
                    continue
                # case 2b
                else: 
                    if alt.startswith("-"):
                        pos_diff_mut.append([s, args.region, pos, q_seq[pos-1], "-"+q_seq[pos: pos+len(alt)-1], "", "", "", "", "", ""])
                    else:
                        if gff_feature:
                            flag = False
                            for k in range(len(gff_feature)):
                                ref_codon, alt_codon, ref_aa, alt_aa, pos_aa = get_codon_aa(q_seq, sample_dic[s], pos, local_start[k], local_end[k], global_start[k], global_end[k]) 
                                # synonymous mutation or unkown aa
                                if ref_aa == alt_aa: 
                                    if flag == False:
                                        pos_diff_mut.append([s, args.region, pos, q_seq[pos-1], alt, "", "", "", "", "", ""])
                                        flag = True
                                    else:
                                        continue
                                else:
                                    pos_diff_mut.append([s, args.region, pos, q_seq[pos-1], alt, gff_feature[k], ref_codon, alt_codon, ref_aa, alt_aa, pos_aa])
                        else:
                            pos_diff_mut.append([s, args.region, pos, q_seq[pos-1], alt, "", "", "", "", "", ""])
            # case 2c
            else: 
                if gff_feature:
                    flag = False
                    for k in range(len(gff_feature)):                        
                        ref_codon, alt_codon, ref_aa, alt_aa, pos_aa = get_codon_aa(q_seq, sample_dic[s], pos, local_start[k], local_end[k], global_start[k], global_end[k]) 
                        # synonymous mutation or unkown aa
                        if ref_aa == alt_aa: 
                            if flag == False:
                                pos_diff_mut.append([s, args.region, pos, q_seq[pos-1], bg_seq[pos-1], "", "", "", "", "", ""])
                                flag = True
                            else:
                                continue
                        else:
                            pos_diff_mut.append([s, args.region, pos, q_seq[pos-1], bg_seq[pos-1], gff_feature[k], ref_codon, alt_codon, ref_aa, alt_aa, pos_aa])
                else:
                    pos_diff_mut.append([s, args.region, pos, q_seq[pos-1], bg_seq[pos-1], "", "", "", "", "", ""])

    pos_diff_mut = pd.DataFrame(pos_diff_mut, columns=header)

    pd.concat([pos_same_mut, pos_diff_mut], axis=0) \
        .sort_values(["sra", "pos"]) \
        .to_csv(args.o, sep="\t", header=header, index=False)


def main():
    parser = argparse.ArgumentParser(description="Extract differences between Hu-1 and other references to infer mutations on other refences.")
    parser.add_argument("-m", help="Mutation file path in TSV format.", required=True)
    parser.add_argument("-d", help="Helper mutation file path in TSV format to deal with leading and trailing dels.", required=True)
    parser.add_argument("-a", help="Alignment file path from gofasta in FASTA format, including bakground genome and query genomes.", required=True)
    parser.add_argument("--gff", help="GFF file to extract gene features and sra.", required=True)
    parser.add_argument("-o", help="Output name for new mutation file in TSV format.", default="mutations.tsv")
    parser.add_argument("--bg", help="ID of background genome.", required=True, type=str)
    parser.add_argument("--query", help="ID of query genom.", required=True, type=str)
    parser.add_argument("--n_ref", help="Number of reference genomes.", required=True, type=int)
    parser.add_argument("--region", help="Region of mutations.", required=True, type=str)
    args = parser.parse_args()

    translate_mutations(args)

    
if __name__ == '__main__':
    main()

# python /home/yutianc/bjorn_rep/bin/translate_mutations.py \
#   -m /home/yutianc/bjorn_rep/output/Hu1/mm/mutations.tsv \
#   -a /home/yutianc/bjorn_rep/output/Hu1/mm/alignment.fasta \
#   --gff /home/yutianc/bjorn_rep/data/Hu1-BA/NC_045512.2.gff \
#   --bg NC_045512.2 \
#   --query NC_045512.2_BA.1 \
#   --region NC_045512.2 \
#   -d /home/yutianc/bjorn_rep/output/Hu1/mm/del_helper.tsv \
#   --n_ref 3

# python /home/yutianc/bjorn_rep/bin/translate_mutations.py \
#   -m /home/yutianc/bjorn_rep/output/PB2/mm/mutations.tsv \
#   -a /home/yutianc/bjorn_rep/output/PB2/mm/alignment.fasta \
#   --gff /home/yutianc/bjorn_rep/data/PB2-DMS/PP755596.1.gff \
#   --bg PP755596.1_cds_XAJ25426.1_1 \
#   --query CY018884.1_cds_ABM21959.1_1 \
#   --region PB2 \
#   -d /home/yutianc/bjorn_rep/output/PB2/mm/del_helper.tsv \
#   --n_ref 2