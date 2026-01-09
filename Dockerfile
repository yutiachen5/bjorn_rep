FROM golang:1.23-alpine AS builder

RUN apk add --no-cache git

WORKDIR /go/src/github.com/virus-evolution/gofasta
RUN git clone https://github.com/rrouz/gofasta.git . && \
    git checkout 0311cec && \
    go build -o /go/bin/gofasta .

FROM condaforge/mambaforge:23.11.0-0

# copy gofasta from go image to the second image
# copy from host machine by default, using --from=builder to copy files from base image
COPY --from=builder /go/bin/gofasta /usr/local/bin/gofasta

RUN apt-get update && apt-get install -y --no-install-recommends \
    grep \
    sed \
    curl \
    bash \
    openjdk-17-jre-headless \
    minimap2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install nf
RUN curl -fsSL https://get.nextflow.io | bash \
    && mv nextflow /usr/local/bin/nextflow \
    && chmod +x /usr/local/bin/nextflow \
    && nextflow self-update 23.04.0 \
    && rm -rf /root/.nextflow

# create conda env
COPY env.yml /tmp/env.yml 

RUN mamba env create -f /tmp/env.yml && \
    mamba clean -a -y && \
    rm /tmp/env.yml

# coopy files to container
# WORKDIR /workspace
# COPY . .

# set env variable for conda
ENV PATH=/opt/conda/envs/bjorn/bin:$PATH

# cmd to run the pipeline
SHELL ["conda", "run", "-n", "bjorn", "/bin/bash", "-c"]


