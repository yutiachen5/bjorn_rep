replication of bjorn pipeline: https://github.com/andersen-lab/bjorn-general/tree/master

current workflow:


1. randomly sample 100 files from consensus sequence (HCoV-19: https://github.com/andersen-lab/HCoV-19-Genomics)
 - some files have the sequence split into multiple lines
 - others keep a single line

2. extract reference genome from NC_045512.2
 - ref fasta downloaded from https://github.com/andersen-lab/bjorn/blob/main/data/reference.fasta (single line)
 - sequence.fasta from Karthik in slack (multiple lines)

3. alignment using minimap2

4. variant calling with gofasta
 - gb file downloaded from: https://github.com/andersen-lab/bjorn/blob/main/data/reference.gb
 - no option for append-codons

5. file format conversion with gofasta
 - tsv file 
 - ref_nuc vs ref_codon?


current output: /output/mutations.tsv

current cmd to run:
```
git clone https://github.com/andersen-lab/HCoV-19-Genomics
cd bjorn_rep
nextflow run main.nf --ref /path/to/ref/fasta --fasta-path /path/to/input/fasta --gb_dir /path/to/GenBank/annotation --metadata /path/to/metadata
```

