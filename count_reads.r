#!/usr/bin/env Rscript
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0-CBI")
args = commandArgs(trailingOnly=T)

demuxed_folder = args[1]

all.files = list.files(paste0(demuxed_folder,"/",list.files(demuxed_folder,pattern="_S"),"/trimmed_demuxed"),pattern="R1",full.names=T)

fastq.summary = ShortRead::countFastq(all.files)

fastq.reads = data.frame(filename = rownames(fastq.summary),
                            reads = fastq.summary$records)

fastq.reads$locus = sapply(strsplit(fastq.reads$filename,"_PARA"),"[",1)
fastq.reads$sample = sapply(strsplit(fastq.reads$filename,"_PARA"),"[",2)


fastq.reads = fastq.reads[,c("sample","locus","reads")]
