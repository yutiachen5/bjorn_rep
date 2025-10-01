# base with go installed
FROM golang:1.23-alpine AS builder

RUN apk add --no-cache git

# set working dir
# WORKDIR /src

# install gofasta
RUN git clone https://github.com/virus-evolution/gofasta.git && \
    cd gofasta && \
    go build -o /usr/local/bin/gofasta .

# base image with conda and mamba
FROM condaforge/mambaforge:23.11.0-0

# copy gofasta from go image to the second image
# copy from host machine by default, using --from=builder to copy files from base image
COPY --from=builder /usr/local/bin/gofasta /usr/local/bin/gofasta


# install jave, using no-install-recommends to avoid beubg asked for the time zone
RUN apt-get update && \
    apt-get install -y \
    openjdk-17-jre-headless \
    curl \
    minimap2 \
    grep \
    mafft \
    --no-install-recommends 

# install nf
RUN curl -fsSL https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow && \
    rm -rf /root/.nextflow

# create conda env
COPY env.yml /tmp/env.yml 
RUN mamba env create -n bjorn -f /tmp/env.yml && \
    mamba clean -a -y && \
    rm /tmp/env.yml

# coopy files to container
WORKDIR /workspace
COPY . .

# set env variable for conda
ENV PATH=/opt/conda/envs/bjorn/bin:$PATH

# cmd to run the pipeline
SHELL ["conda", "run", "-n", "env", "/bin/bash", "-c"]

CMD ["nextflow", "run", "main.nf", \
    "--ref_file", "/workspace/data/NC_045512.2.fasta", \
    "--fasta_dir", "/workspace/data/consensus_sequences", \
    "--gff_file", "/workspace/data/NC_045512.2.gff", \
    "--query", "/workspace/data/BA.1_and_BA.2.fa", \
    "-c", "nf.config"]

