replication of bjorn pipeline: https://github.com/andersen-lab/bjorn-general/tree/master

current workflow:

1. extract data from consensus sequence (HCoV-19: https://github.com/andersen-lab/HCoV-19-Genomics)

2. alignment using minimap2
   - reference:NC_045512.2
    - Hu-1
    - Ba-1

3. mutation calling 
   - gofasta sam variants
   - manual calling

4. Translate mutations from current reference genome to other reference genomes
   - [BA.1 mutations](output/NC_045512.2_BA.1_mutations.tsv)
   - [BA.2 mutations](output/NC_045512.2_BA.2_mutations.tsv)
   - Limitations:
     - only works when len(bg.seq) == len(q.seq)


current cmd to run with docker:
```
docker build -t bjorn .
docker run bjorn
```


with manual mutation calling
```
nextflow run main.nf \
   --fasta_dir $PWD/data/Hu1-BA/consensus_sequences \
   --ref_file $PWD/data/Hu1-BA/ref_all.fasta \
   --gff_file $PWD/data/Hu1-BA/NC_045512.2.gff \
   --translate_mutations true \
   --ref_id NC_045512.2 \
   --query_id NC_045512.2_BA.1,NC_045512.2_BA.2 \
   --region NC_045512.2 \
   --nsamples 1000 \
   --outdir $PWD/output/Hu1/ \
   -c nf.config
```

with gofasta mutation calling (no translation)
```
nextflow run main.nf \
   --gofasta true \
   --fasta_dir $PWD/data/Hu1-BA/consensus_sequences \
   --ref_file $PWD/data/Hu1-BA/NC_045512.2.fasta \
   --ref_id NC_045512.2 \
   --gff_file $PWD/data/Hu1-BA/NC_045512.2.gff \
   --translate_mutations false \
   --nsamples 1000 \
   --region NC_045512.2 \
   --outdir $PWD/output/Hu1/gofasta \
   -c nf.config
```

on PB-2 data:
```
nextflow run main.nf \
   --fasta_dir $PWD/data/PB2-DMS/PB2_samples/ \
   --ref_file $PWD/data/PB2-DMS/PP755596.1.fasta \
   --gff_file $PWD/data/PB2-DMS/PP755596.1.gff \
   --query_file $PWD/data/PB2-DMS/CY018884.1.fasta \
   --region PP755596.1 \
   --outdir $PWD/output/PB2 \
   --ref_id PP755596.1_cds_XAJ25426.1_1 \
   --query_id CY018884.1_cds_ABM21959.1_1 \
   --sampling false \
   -c nf.config
```

