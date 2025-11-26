#!/usr/bin/env python

import re
import argparse
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO

def parse_alignment_file(args): 
    records = list(SeqIO.parse(args.a, "fasta"))

    seq_len = set([len(record.seq) for record in records])
    assert len(seq_len) == 1, "Sequence length mismatch!"
    length = seq_len.pop()

    ref_records = records[:args.n_ref] # refs are put at the top of alignment file
    sample_records = records[args.n_ref: ]

    ref_dic = {record.id: record.seq for record in ref_records if record.id in [args.bg, args.query]}
    sample_dic = {record.id: record.seq for record in sample_records}

    references = []
    for i in range(length): 
        if ref_dic[args.bg][i] != ref_dic[args.query][i]:
            references.append({"pos": i+1, "bg.ref":ref_dic[args.bg][i], "q.ref":ref_dic[args.query][i]}) # 1-based index

    return pd.DataFrame(references), ref_dic, sample_dic

def get_gff_feature(gff, pos):
    # GFF intervals are 1-based, inclusive
    # pos is 1-based index
    match = gff[(gff["start"] <= pos) & (pos <= gff["end"])]

    if len(match) > 0:
        gff_feature = match["GFF_FEATURE"].tolist()
        start = match["start"].tolist()
        end = match["end"].tolist()
    else:
        gff_feature, start, end = [], [], []

    return gff_feature, start, end

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
            if cols[2] == 'CDS': # cds only
                region = cols[0]
                start, end = int(cols[3]), int(cols[4])
                gff_feature = re.search(r"(?<=cds-)[^;\s]+", cols[-1]).group()
                gff.append({"region":region, "start": start, "end": end, "GFF_FEATURE": gff_feature})

    return pd.DataFrame(gff)

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

    pos_aa = ((pos - start) // 3) + 1

    return ref_codon, alt_codon, ref_aa, alt_aa, pos_aa

def translate_mutations(args):
    mutation = pd.read_csv(args.m, sep="\t", header=0)
    helper = pd.read_csv(args.d, sep="\t", header=0)
    deletion = extract_del(mutation, helper)

    ref_variants, ref_dic, sample_dic = parse_alignment_file(args) # ref_variants: {"pos": 2, "Hu1": "A", "BA1": "T"}
    gff = parse_gff(args.gff)
    bg_seq, q_seq = ref_dic[args.bg], ref_dic[args.query]

    header = ["sra", "region", "pos", "ref", "alt", "GFF_FEATURE", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "pos_aa"]

    # case1: bg.ref == q.ref
        # case1a: mut_type in (INS, SNP), appended
        # case1b: mut_type == DEL, update ref seq
    # case2: bg.ref != q.ref
        # case2a: bg.ref != alt == q.ref, skip
        # case2b: bg.ref != alt != q.ref (alt can be del), update ref = q.ref, alt = alt
        # case2c: bg.ref == alt != q.ref, add ref = q.ref, alt = bg.ref

    pos_diff = ref_variants["pos"].to_list()

    pos_same = mutation[~mutation["pos"].isin(pos_diff)]
    # case1a
    pos_same_ins_snp = pos_same[~pos_same["alt"].str.startswith("-")]
    # case1b
    pos_same_del = pos_same[pos_same["alt"].str.startswith("-")].copy().reset_index(drop=True)
    if not pos_same_del.empty:
        pos_same_del["alt"] = pos_same_del.apply(
            lambda row: "-" + ''.join(
                q_seq[row["pos"]: row["pos"] + len(row["alt"]) - 1]
            ),
            axis=1,
        )

    pos_same_mut = pd.concat([pos_same_ins_snp, pos_same_del], axis=0)


    pos_diff_mut = []

    for i in range(len(pos_diff)):
        pos = pos_diff[i] # 1-based index
        gff_feature, start, end = get_gff_feature(gff, pos) 
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
                        pos_diff_mut.append([s, args.region, str(pos), q_seq[pos-1], "-"+q_seq[pos: pos+len(alt)-1], "", "", "", "", "", ""])
                    else:
                        for k in range(len(gff_feature)):
                            ref_codon, alt_codon, ref_aa, alt_aa, pos_aa = get_codon_aa(q_seq, sample_dic[s], pos, start[k], end[k]) 
                            if ref_aa == alt_aa: # synonymous mutation
                                pos_diff_mut.append([s, args.region, str(pos), q_seq[pos-1], alt, "", "", "", "", "", ""])
                            else:
                                pos_diff_mut.append([s, args.region, str(pos), q_seq[pos-1], alt, gff_feature[k], ref_codon, alt_codon, ref_aa, alt_aa, pos_aa])
            # case 2c
            else: 
                for k in range(len(gff_feature)):                        
                    ref_codon, alt_codon, ref_aa, alt_aa, pos_aa = get_codon_aa(q_seq, sample_dic[s], pos, start[k], end[k]) 
                    if ref_aa == alt_aa: # synonymous mutation
                        pos_diff_mut.append([s, args.region, str(pos), q_seq[pos-1], bg_seq[pos-1], "", "", "", "", "", ""])
                    else:
                        pos_diff_mut.append([s, args.region, str(pos), q_seq[pos-1], bg_seq[pos-1], gff_feature[k], ref_codon, alt_codon, ref_aa, alt_aa, pos_aa])

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

# python /home/eleanor124/projects/bjorn_rep/bin/translate_mutations.py \
#   -m /home/eleanor124/projects/bjorn_rep/output/Hu1/mm/mutations.tsv \
#   -a /home/eleanor124/projects/bjorn_rep/output/Hu1/mm/alignment.fasta \
#   --gff /home/eleanor124/projects/bjorn_rep/data/Hu1-BA/NC_045512.2.gff \
#   --bg NC_045512.2 \
#   --query NC_045512.2_BA.1 \
#   --region NC_045512.2 \
#   -d /home/eleanor124/projects/bjorn_rep/output/Hu1/mm/del_helper.tsv \
#   --n_ref 3