#!/usr/bin/env bash
set -e

module load openjdk

currentDate="$(date +%m%d%Y_%H%M%S)"

mkdir /wynton/scratch/aarandad/mad4hatter/results/$1_$currentDate
cp nextflow.config /wynton/scratch/aarandad/mad4hatter/results/$1_$currentDate/nextflow.config



git log -1 --format="%H" > /wynton/scratch/aarandad/mad4hatter/results/$1_$currentDate/git.txt
git rev-parse --abbrev-ref HEAD >> /wynton/scratch/aarandad/mad4hatter/results/$1_$currentDate/git.txt

echo $1_$currentDate
NXF_VER=22.11.0-edge nextflow -log /wynton/scratch/aarandad/mad4hatter/log/$1_$currentDate \
	run main.nf -c conf/custom.config -profile sge,apptainer \
	--sequencer nextseq \
	--refseq_fasta /wynton/home/eppicenter/aarandad/mad4hatter_new/mad4hatter/resources/v4/v4_reference.fasta \
	--target v4 \
	--readDIR /wynton/scratch/aarandad/mad4hatter/data/$1 \
	--outDIR /wynton/scratch/aarandad/mad4hatter/results/$1_$currentDate \
	-w /wynton/scratch/aarandad/mad4hatter/work/$1_$currentDate 
