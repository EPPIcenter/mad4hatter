#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)

demuxed_folder = args[1]

fastq.summary = ShortRead::countFastq(demuxed_folder,pattern = "R1")

fastq.reads = data.frame(filename = rownames(fastq.summary),
                            reads = fastq.summary$records)

fastq.reads$locus = sapply(strsplit(fastq.reads$filename,"_trimmed"),"[",1)

fastq.reads = fastq.reads[,c("locus","reads")]

write.table(fastq.reads, file=args[2], quote=FALSE, sep='\t')

