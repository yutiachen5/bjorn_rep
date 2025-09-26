#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SRATOOLS_PREFETCH {
    input:
    tuple val(meta), val(id)

    output:
    tuple val(meta), path(id), emit: sra

    shell:
    args = '5 1 100'
    template('retry_with_backoff.sh')

}