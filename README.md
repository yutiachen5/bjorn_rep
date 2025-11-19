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
   --translate_mutations true \
   --nsamples 1000 \
   --fasta_dir $PWD/data/Hu1-BA/consensus_sequences \
   --ref_file $PWD/data/Hu1-BA/NC_045512.2.fasta \
   --query_ref_file $PWD/data/Hu1-BA/BA.1_and_BA.2.fa \
   --gff_file $PWD/data/Hu1-BA/NC_045512.2.gff \
   --region NC_045512.2 \
   --outdir $PWD/output/Hu1/mm \
   --chunk_size 100 \
   -c nf.config
```

with gofasta mutation calling (no translation) - somthing weird in BA1 is happening
```
nextflow run main.nf \
   --gofasta true \
   --translate_mutations false \
   --nsamples 1000 \
   --fasta_dir $PWD/data/Hu1-BA/consensus_sequences \
   --ref_file $PWD/data/Hu1-BA/NC_045512.2.fasta \
   --gff_file $PWD/data/Hu1-BA/NC_045512.2.gff \
   --region NC_045512.2 \
   --outdir $PWD/output/Hu1/gf_Hu1 \
   -c nf.config
```

on PB-2 - manual mutation calling:
```
nextflow run main.nf \
   --sampling false \
   --fasta_dir $PWD/data/PB2-DMS/PB2_samples/ \
   --ref_file $PWD/data/PB2-DMS/PP755596.1.fasta \
   --query_ref_file $PWD/data/PB2-DMS/CY018884.1.fasta \
   --gff_file $PWD/data/PB2-DMS/PP755596.1.gff \
   --region PB2 \
   --outdir $PWD/output/PB2/mm \
   -c nf.config
```

on PB-2 - gofasta variants:
```
nextflow run main.nf \
   --gofasta true \
   --translate_mutations false \
   --sampling false \
   --fasta_dir $PWD/data/PB2-DMS/PB2_samples/ \
   --ref_file $PWD/data/PB2-DMS/PP755596.1.fasta \
   --gff_file $PWD/data/PB2-DMS/PP755596.1.gff \
   --region PB2 \
   --outdir $PWD/output/PB2/gf_PP \
   -c nf.config
```
