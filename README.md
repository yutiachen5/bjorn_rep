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

6. Translate mutations from current reference genome to other reference genomes
   - file type conversion: sam to fasta using gofasta tomultialigned
   - [BA.1 mutations](output/NC_045512.2_BA.1_mutations.tsv)
   - [BA.2 mutations](output/NC_045512.2_BA.2_mutations.tsv)
   - Limitations:
     - only works when len(bg.seq) == len(q.seq)
     - ignoring insertions and deletions

7. Run pipeline on on BA samples (38901 in total)

current cmd to run with docker:
```
docker build -t bjorn .
docker run bjorn
```

for all BA samples:
```
nextflow run main.nf \
   --fasta_dir $PWD/data/consensus_sequences \
   --ref_file $PWD/data/NC_045512.2.fasta \
   --gff_file $PWD/data/NC_045512.2.gff \
   --query $PWD/data/BA.1_and_BA.2.fa \
   --region NC_045512.2 \
   --lineage_file $PWD/data/lineage_report.csv \
   --l BA \
   -c nf.config 

```

100 samples from random selection:
```
nextflow run main.nf \
   --fasta_dir $PWD/data/consensus_sequences \
   --ref_file $PWD/data/NC_045512.2.fasta \
   --gff_file $PWD/data/NC_045512.2.gff \
   --query $PWD/data/BA.1_and_BA.2.fa \
   --region NC_045512.2 \
   -c nf.config 
```

