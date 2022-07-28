.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0-CBI")

library(tidyverse)
library(dada2)

args = commandArgs(trailingOnly=T)

#treat.no.overlap.differently = T
#trimmed_path = "~/Data/test_times/results/cutadapt_TES_quarter_demux_new"
#ampliconFILE = "v3_amplicon_info.tsv"

numargs=length(args)

trimmed_path=args[c(1:(numargs-3))]
ampliconFILE=args[(numargs-2)]
treat.no.overlap.differently=args[(numargs-1)]
dada2RData=args[numargs]

print(trimmed_path)

fnFs <- sort(list.files(path=trimmed_path, pattern="_R1.fastq.gz", recursive=T, full.names = TRUE))
fnRs <- sort(list.files(path=trimmed_path, pattern="_R2.fastq.gz", recursive=T, full.names = TRUE))

#fnFs = sort(grep("_R1.fastq.gz",trimmed_path,value=T))
#fnRs = sort(grep("_R2.fastq.gz",trimmed_path,value=T))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
print(sample.names)

filtFs <- paste0(trimmed_path,"/filtered/",sample.names,"_F_filt.fastq.gz")
filtRs <- paste0(trimmed_path,"/filtered/",sample.names,"_R_filt.fastq.gz")

#filtFs <- paste0(dirname(trimmed_path),"/filtered/",sample.names,"_F_filt.fastq.gz")
#filtRs <- paste0(dirname(trimmed_path),"/filtered/",sample.names,"_R_filt.fastq.gz")


names(filtFs) <- sample.names
names(filtRs) <- sample.names
# 

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

if(treat.no.overlap.differently == T){
  
  # extract name of locus from filename
  amplicon.info = data.frame(
    names = sapply(strsplit(names(dadaFs),"_trimmed"),"[",1),
    sum.mean.length.reads  = sapply(sapply(unlist(lapply(dadaFs,"[","sequence"),recursive = F),nchar),mean)+
      sapply(sapply(unlist(lapply(dadaRs,"[","sequence"),recursive = F),nchar),mean)
  ) %>% 
    left_join(read.table(ampliconFILE,header=T) %>% 
                select(amplicon,ampInsert_length),
              by=c("names"="amplicon"))
  rownames(amplicon.info)=names(dadaFs)
  
#  sample.names.no.overlap = sample.names[as.numeric(sapply(strsplit(sample.names ,"-"),"[",3)) - as.numeric(sapply(strsplit(sample.names ,"-"),"[",2))>275]
#  sample.names.overlap = sample.names[as.numeric(sapply(strsplit(sample.names ,"-"),"[",3)) - as.numeric(sapply(strsplit(sample.names ,"-"),"[",2))<276]
  
  mergers.overlap = mergePairs(dadaFs[names(dadaFs) %in% rownames(amplicon.info %>% filter(ampInsert_length<sum.mean.length.reads))],
                               derepFs[names(dadaFs) %in% rownames(amplicon.info %>% filter(ampInsert_length<sum.mean.length.reads))], 
                               dadaRs[names(dadaFs) %in% rownames(amplicon.info %>% filter(ampInsert_length<sum.mean.length.reads))], 
                               derepRs[names(dadaFs) %in% rownames(amplicon.info %>% filter(ampInsert_length<sum.mean.length.reads))], 
                               verbose=TRUE, 
                               justConcatenate=FALSE, 
                               trimOverhang = TRUE,
                               minOverlap=10,
                               maxMismatch=1)
  
  mergers.no.overlap = mergePairs(dadaFs[names(dadaFs) %in% rownames(amplicon.info %>% filter(ampInsert_length>=sum.mean.length.reads))],
                                  derepFs[names(dadaFs) %in% rownames(amplicon.info %>% filter(ampInsert_length>=sum.mean.length.reads))], 
                                  dadaRs[names(dadaFs) %in% rownames(amplicon.info %>% filter(ampInsert_length>=sum.mean.length.reads))], 
                                  derepRs[names(dadaFs) %in% rownames(amplicon.info %>% filter(ampInsert_length>=sum.mean.length.reads))], 
                                  verbose=TRUE, 
                                  justConcatenate=TRUE)
  
  mergers = c(mergers.overlap,mergers.no.overlap)[names(dadaFs)]
  }else{
    mergers=  mergePairs(dadaFs,
                         derepFs, 
                         dadaRs, 
                         derepRs, 
                         verbose=TRUE, 
                         justConcatenate=FALSE, 
                         trimOverhang = TRUE,
                         minOverlap=10,
                         maxMismatch=1)
  }

seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

save(out,dadaFs,dadaRs,mergers,seqtab,seqtab.nochim, file = dada2RData)
