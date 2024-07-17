#!/usr/bin/env bash

DATADIR="/data/sequencing/nextseq-investigation"

while read line; do
	echo "Processing ${line}..."
	nextflow -q run main.nf --readDIR ${DATADIR}/${line} -profile docker --target v4 -c conf/custom.config --QC_only -w "nextseq-inv/work_${line}" --outDIR "nextseq-inv/results_${line}"
done < <(ls -1 $DATADIR)