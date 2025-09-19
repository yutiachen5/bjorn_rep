replication of bjorn pipeline: https://github.com/andersen-lab/bjorn-general/tree/master

current workflow:


1. randomly sample 100 files from consensus sequence (HCoV-19: https://github.com/andersen-lab/HCoV-19-Genomics)

2. extract reference genome from NC_045512.2

3. alignment using minimap2

4. variant calling with gofasta

5. file format conversion with gofasta


output: tsv file

