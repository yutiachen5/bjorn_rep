#!/usr/bin/env python

import os
import random
import argparse


def sampling(fasta_dir, n):
    all_samples = [os.path.join(fasta_dir, f) \
                    for f in os.listdir(fasta_dir) if os.path.isfile(os.path.join(fasta_dir, f))]
    selected_fasta = random.sample(all_samples, min(n, len(all_samples)))

    with open("selected_fasta.txt", "w") as out:
        for f in selected_fasta:
            out.write(f + "\n")
    
def main():
    parser = argparse.ArgumentParser(description="Randomly select n samples from all input fasta files.")
    parser.add_argument("--fasta_dir", help="Input fasta dir", required=True)
    # parser.add_argument("--output", default="mutations.tsv", help="Output TSV file name.")
    parser.add_argument("-n", default=100, help="Number of samples.")

    args = parser.parse_args()

    selected_fasta_ls = sampling(args.fasta_dir, args.n)

if __name__ == "__main__":
    main() 

