# Bjorn

This pipeline is the rebuilt of [Bjorn-general](https://github.com/andersen-lab/bjorn-general.git). The Bjorn-general pipeline is a Nextflow pipeline designed for analysis of H5N1 influenza sequence data, integrating deep mutational scanning (DMS), and obtaining mutations at both a host and population level.

## Key Features

*   **Variant Calling:** Identifies sequence variants relative to a reference genome and translates mutations to other reference geomes in DMS.

## Data Setup

*   **Reference Genome:** A default `reference.fasta` file is included under [SC2_Reference](data/sc2)
*   **Sequence Data:** The pipeline was only tested on all San Diego samples from [andersen-lab/avian-influenza](https://github.com/andersen-lab/HCoV-19-Genomics).



current cmd to run with docker:
```
docker build -t bjorn .
docker run bjorn
```

current cmd to run in terminal:

manual mutation calling
```
nextflow run main.nf \
   --translate_mutations true \
   --nsamples 1000 \
   --fasta_dir $PWD/data/sc2/consensus_sequences \
   --ref_file $PWD/data/sc2/NC_045512.2.fasta \
   --query_ref_file $PWD/data/sc2/escape_all.fa \
   --gff_file $PWD/data/sc2/NC_045512.2.gff \
   --region NC_045512.2 \
   --outdir $PWD/output/SC2/escape_all \
   --chunk_size 20 \
   -c nf.config
```

gofasta mutation calling (no translation) 
```
nextflow run main.nf \
   --gofasta true \
   --translate_mutations false \
   --nsamples 1000 \
   --fasta_dir $PWD/data/sc2/consensus_sequences \
   --ref_file $PWD/data/sc2/NC_045512.2.fasta \
   --gff_file $PWD/data/sc2/NC_045512.2.gff \
   --region NC_045512.2 \
   --outdir $PWD/output/SC2/gf_Hu1 \
   --chunk_size 20 \
   -c nf.config
```

PB-2 - manual mutation calling:
```
nextflow run main.nf \
   --sampling false \
   --translate_mutations true \
   --fasta_dir $PWD/data/pb2/PB2_samples/ \
   --ref_file $PWD/data/pb2/PP755596.1.fasta \
   --query_ref_file $PWD/data/pb2/CY018884.1.fasta \
   --gff_file $PWD/data/pb2/PP755596.1.gff \
   --region PB2 \
   --outdir $PWD/output/PB2/mm \
   --chunk_size 100 \
   -c nf.config
```

PB-2 - gofasta variants:
```
nextflow run main.nf \
   --gofasta true \
   --translate_mutations false \
   --sampling false \
   --fasta_dir $PWD/data/pb2/PB2_samples/ \
   --ref_file $PWD/data/pb2/PP755596.1.fasta \
   --gff_file $PWD/data/pb2/PP755596.1.gff \
   --region PB2 \
   --outdir $PWD/output/PB2/gf_PP \
   --chunk_size 100 \
   -c nf.config
```

Docker:

```
docker run --rm -it \
  -v "$PWD":/workspace -w /workspace \
  -e NXF_ANSI_LOG=true \
  -e TERM=xterm-256color \
  bjorn \
  bash -lc 'conda run -n bjorn nextflow run main.nf -profile docker --translate_mutations true --sampling false --fasta_dir /workspace/data/sc2/consensus_sequences --ref_file /workspace/data/sc2/NC_045512.2.fasta --query_ref_file /workspace/data/sc2/escape_all.fa --gff_file /workspace/data/sc2/NC_045512.2.gff --region NC_045512.2 --outdir /workspace/mutations/SC2/escape_all/ --chunk_size 20 -c nf.config'
```