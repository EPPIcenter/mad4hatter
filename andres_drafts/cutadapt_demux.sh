#!/bin/bash

mkdir results/cutadapt_TES_quarter_demux_new

cutadapt \
    --action=trim \
    --discard-untrimmed \
    -g file:"/home/andres/Data/amplicon_sequencing/fastas/v3_fwd.fasta" \
    -G file:"/home/andres/Data/amplicon_sequencing/fastas/v3_rev.fasta" \
    --pair-adapters \
    -e 0 \
    --no-indels \
    -o results/cutadapt_TES_quarter_demux_new/{name}_trimmed_R1.fastq.gz \
    -p results/cutadapt_TES_quarter_demux_new/{name}_trimmed_R2.fastq.gz \
    data/raw_fastq/TES_quarter_R1_001.fastq.gz \
    data/raw_fastq/TES_quarter_R2_001.fastq.gz 
