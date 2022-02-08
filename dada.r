#!/usr/bin/env Rscript

args = commandArgs(trainingOnly=T)

library(dada2)

trimmed_path = args[1]
fnFs <- sort(list.files(path=trimmed_path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path=trimmed_path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

filtFs <- paste0(trimmed_path,"/filtered/",sample.names,"_F_filt.fastq.gz")
filtRs <- paste0(trimmed_path,"/filtered/",sample.names,"_R_filt.fastq.gz")

names(filtFs) <- sample.names
names(filtRs) <- sample.names


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
              maxN=0, maxEE=c(2,2), truncQ=c(5,5), rm.phix=TRUE,
              compress=TRUE, multithread=TRUE,
              trimRight = c(0,0),trimLeft = 1, minLen=75,matchIDs=TRUE) 
head(out)


errF <- learnErrors(filtFs[out[,2]>0], multithread=TRUE,MAX_CONSIST=10,randomize=TRUE)
errR <- learnErrors(filtRs[out[,2]>0], multithread=TRUE,MAX_CONSIST=10,randomize=TRUE)

derepFs <- derepFastq(filtFs[out[,2]>0], verbose = TRUE)
derepRs <- derepFastq(filtRs[out[,2]>0], verbose = TRUE)

dadaFs <- dada(derepFs, err=errF, selfConsist=TRUE, multithread=TRUE, verbose=TRUE, OMEGA_A=1e-120)
dadaRs <- dada(derepRs, err=errR, selfConsist=TRUE, multithread=TRUE, verbose=TRUE, OMEGA_A=1e-120)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, justConcatenate=FALSE, trimOverhang = TRUE,minOverlap=10,maxMismatch=1)

seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

save(out,dadaFs,dadaRs,mergers,seqtab,seqtab.nochim, file = args[2])
