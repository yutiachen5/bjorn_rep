#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SRATOOLS_PREFETCH {
    input:
    tuple val(meta), val(id)

    output:
    path 


}