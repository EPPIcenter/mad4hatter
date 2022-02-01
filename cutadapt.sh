#!/bin/bash

mkdir results/cutadapt_TES

cutadapt \
    --action=trim \
    --discard-untrimmed \
    -g file:"data/fastas/v3_fwd.fasta" \
    -G file:"data/fastas/v3_rev.fasta" \
    --pair-adapters \
    -e 0 \
    --no-indels \
    -o results/cutadapt_TES/trimmed_R1.fastq.gz \
    -p results/cutadapt_TES/trimmed_R2.fastq.gz \
    data/raw_fastq/TES_R1_001.fastq.gz \
    data/raw_fastq/TES_R2_001.fastq.gz 
