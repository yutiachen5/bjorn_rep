#!/usr/bin/env python

import os
import pandas as pd
import random
import argparse

random.seed(42)

def sampling(args):
    # select samples from certain lineage
    if getattr(args, "l", None):
        df = pd.read_csv(args.lineage_file)
        df["prefix"] = df['lineage'].str.extract("([A-Za-z]+)(?=\.\d+)")
        df["prefix"] = df["prefix"].fillna(df["lineage"])
        selected_fasta = list(df.loc[df["prefix"] == args.l, :]["taxon"])

        with open("fasta_files.txt", "w") as outfile:
            # outfile.write(args.ref + "\n")
            for f in selected_fasta:
                if os.path.isfile(os.path.join(args.fasta_dir, f+".fasta")):
                    path = os.path.join(args.fasta_dir, f+".fasta") 
                    outfile.write(path + "\n")
    # random sampling
    else:
        all_samples = [os.path.join(args.fasta_dir, f) \
                        for f in os.listdir(args.fasta_dir) if os.path.isfile(os.path.join(args.fasta_dir, f))]
        selected_fasta = random.sample(all_samples, min(args.n, len(all_samples)))

        with open("fasta_files.txt", "w") as outfile:
            # outfile.write(args.ref + "\n")
            for f in selected_fasta:
                outfile.write(f + "\n")
    
def main():
    parser = argparse.ArgumentParser(description="Randomly select n samples from all input fasta files.")
    parser.add_argument("--fasta_dir", help="Input fasta dir", required=True)
    parser.add_argument("-n", help="Number of samples.", default=100, type=int)
    parser.add_argument("-l", help="Lineage name for extracting mutations.")
    parser.add_argument("--lineage_file", help="Lineage file in CVS format.")
    # parser.add_argument("--ref", help="Path of reference genome in FASTA format")

    args = parser.parse_args()

    sampling(args)

if __name__ == "__main__":
    main() 

# python /home/eleanor124/projects/bjorn_rep/bin/sampling.py \
#   --lineage_file /home/eleanor124/projects/HCoV-19-Genomics/lineage_report.csv \
#   -l BA \
#   --fasta_dir /home/eleanor124/projects/HCoV-19-Genomics/consensus_sequences/ \
#   --ref /home/eleanor124/projects/bjorn_rep/data/NC_045512.2.fasta
