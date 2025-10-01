replication of bjorn pipeline: https://github.com/andersen-lab/bjorn-general/tree/master

current workflow:

1. randomly sample 100 files from consensus sequence (HCoV-19: https://github.com/andersen-lab/HCoV-19-Genomics)

2. alignment using minimap2
 - reference:NC_045512.2
  - Hu-1
  - Ba-1

4. variant calling with gofasta
 - no option for append-codons

5. file format conversion with gofasta
 - tsv file 


current cmd to run with docker:
```
docker build -t bjorn .
docker run bjorn
```

cmd to run locally:
```
nextflow run main.nf --fasta_dir $PWD/data/consensus_sequences --ref_file $PWD/data/NC_045512.2.fasta --gff_file $PWD/data/NC_045512.2.gff --query /home/eleanor124/projects/bjorn_rep/data/BA.1_and_BA.2.fa -c nf.config
```
