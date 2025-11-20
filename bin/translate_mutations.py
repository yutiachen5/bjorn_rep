#!/usr/bin/env python

import re
import argparse
import pandas as pd
from Bio import SeqIO

def extract_ref_variants(args): 
    records = list(SeqIO.parse(args.a, "fasta"))

    seq_len = set([len(record.seq) for record in records])
    assert len(seq_len) == 1, "Sequence length mismatch!"
    length = seq_len.pop()

    ref_records = records[:args.n_ref] # refs are put at the top of alignment file
    ref_dic = {record.id: record.seq for record in ref_records if record.id in [args.bg, args.query]}

    sample_records = records[args.n_ref:]
    sra = [s.id for s in sample_records]

    references = []
    for i in range(length): 
        if ref_dic[args.bg][i] != ref_dic[args.query][i]:
            references.append({"pos": i+1, "bg.ref":ref_dic[args.bg][i], "q.ref":ref_dic[args.query][i]}) # 1-based index

    return pd.DataFrame(references), ref_dic, sra

def get_gff_feature(pos, gff):
    match = gff[(gff["start"] <= pos) & (pos <= gff["end"])]

    return match["GFF_FEATURE"].tolist() if not match.empty else []

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

    # region covered by leading or trailing del
    helper["start"] = helper["pos"] + 1
    helper["end"] = helper["pos"] + helper["alt"].apply(lambda x: len(x)-1)

    deletion = pd.concat([deletion, helper], axis=0)

    return deletion

def match_del(pos, s, deletion):
    match = deletion[(deletion["start"] <= pos) & (pos <= deletion["end"]) & (deletion["sra"] == s)]

    return True if not match.empty else False

def translate_mutations(args):
    mutation = pd.read_csv(args.m, sep="\t", header=0)
    helper = pd.read_csv(args.d, sep="\t", header=0)
    deletion = extract_del(mutation, helper)

    ref_variants, ref_dic, sample_id = extract_ref_variants(args) # ref_variants: {"pos": 2, "Hu1": "A", "BA1": "T"}
    gff = parse_gff(args.gff)

    bg_seq, q_seq = ref_dic[args.bg], ref_dic[args.query]
    length = len(bg_seq)

    header = ["sra", "region", "pos", "ref", "alt", "GFF_FEATURE", "ref_codon", "alt_codon", "ref_aa", "alt_aa", "pos_aa"]

    # case1: both & old_alt == q.ref, continue
    # case2: both & old_alt != q.ref, write: ref = q.ref, alt = old_ale
    # case3: mut on new ref only (old_alt == bg.ref != q.ref), write for all samples: ref = q.ref, alt = bg.ref 
    # case4: bg.ref == q.ref, write the current mut: ref = old_ref, alt = old_alt

    with open(args.o, "w") as outfile:
        outfile.write("\t".join(header) + "\n")

        for i in range(length):
            if (i+1) not in ref_variants["pos"].to_list(): # case 4: bg.ref == q.ref
                tmp_ins_snp = mutation[(mutation["pos"] == i+1) & (~mutation["alt"].str.startswith("-"))]
                tmp_ins_snp.to_csv(outfile, sep="\t", header=False, index=False)

                tmp_del = mutation[(mutation["pos"] == i+1) & (mutation["alt"].str.startswith("-"))]
                if not tmp_del.empty:
                    tmp_del = tmp_del.copy()   
                    tmp_del["alt"] = tmp_del.apply(
                        lambda row: "-" + q_seq[row["pos"] : row["pos"]+len(row["alt"])-1],
                        axis=1
                    )
                    tmp_del.to_csv(outfile, sep="\t", header=False, index=False)
            else:
                gff_feature = get_gff_feature(i+1, gff) 

                tmp = mutation[mutation["pos"] == i+1]
                tmp_sample = list(tmp["sra"])

                for s in sample_id:
                    del_match = match_del(i+1, s, deletion)
                    if del_match:
                        continue
                    if s in tmp_sample: # alt != bg.ref
                        alt = tmp.loc[tmp["sra"] == s, "alt"].iloc[0]
                        if alt == q_seq[i]: # case 1: alt != bg.ref & alt == q.ref
                            continue
                        else: # case 2: alt != bg.ref & alt != q.ref
                            for g in gff_feature:
                                outfile.write("\t".join([s, args.region, str(i+1), q_seq[i], alt, g]) + "\n")
                    else: # case 3: alt == bg.ref
                        for g in gff_feature:
                            outfile.write("\t".join([s, args.region, str(i+1), q_seq[i], bg_seq[i], g]) + "\n")

            

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