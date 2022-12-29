library(tidyverse)
library(dada2)
library(argparse)

parser <- ArgumentParser(prog="DADA2 Workflow", description='DADA2 workflow script')
parser$add_argument('--trimmed-path', type="character",
                   help='homopolymer threshold to begin masking', nargs='+', required = TRUE)
parser$add_argument('--ampliconFILE', type="character", required = TRUE)
parser$add_argument('--dada2-rdata-output', type="character", required = TRUE)
parser$add_argument('--pool', type="character", default="false")
parser$add_argument('--band-size', type='integer', default=16)
parser$add_argument('--omega-a', type='double', default=1e-120)
parser$add_argument('--concat-non-overlaps', action='store_true')

args <- parser$parse_args()
print(args)

fnFs <- sort(list.files(path=args$trimmed_path, pattern="_R1.fastq.gz", recursive=T, full.names = TRUE))
fnRs <- sort(list.files(path=args$trimmed_path, pattern="_R2.fastq.gz", recursive=T, full.names = TRUE))

#fnFs = sort(grep("_R1.fastq.gz",trimmed_path,value=T))
#fnRs = sort(grep("_R2.fastq.gz",trimmed_path,value=T))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
print(sample.names)

filtFs <- paste0(args$trimmed_path,"/filtered/",sample.names,"_F_filt.fastq.gz")
filtRs <- paste0(args$trimmed_path,"/filtered/",sample.names,"_R_filt.fastq.gz")

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


pool=switch(args$pool,
  "true" = TRUE,
  "false" = FALSE,
  "pseudo" = "pseudo")

dadaFs <- dada(derepFs, err=errF, selfConsist=TRUE, multithread=TRUE, verbose=TRUE, pool=pool, BAND_SIZE=args$band_size, OMEGA_A=args$omega_a)
dadaRs <- dada(derepRs, err=errR, selfConsist=TRUE, multithread=TRUE, verbose=TRUE, pool=pool, BAND_SIZE=args$band_size, OMEGA_A=args$omega_a)

if(args$concat_non_overlaps){
  

  amplicon.info = data.frame(
    names = sapply(strsplit(names(dadaFs),"_trimmed"),"[",1)
  )

  amplicon.info <- amplicon.info %>%
    mutate(names=sapply(str_split(names,'_S(\\d+)'),head,1)) %>% 
    mutate(amplicon=unlist(lapply(str_split(names,'_'), function(x) { paste(x[1:3], collapse = "_") }))) %>%
    inner_join(read.table("v2_amplicon_info.tsv",header=T) %>% 
             select(amplicon, ampInsert_length), by = c("amplicon")) %>%
    # mutate(Pool=sapply(str_split(Amplicon,'-'),tail,1)) %>% 
    # mutate(Amplicon=unlist(lapply(str_split(Amplicon,'-'), function(x) { paste(x[1:2], collapse = "-") })))

    select(amplicon, ampInsert_length) %>%
    mutate(
      sum.mean.length.reads = sapply(sapply(unlist(lapply(dadaFs,"[","sequence"),recursive = F),nchar),mean)+
        sapply(sapply(unlist(lapply(dadaRs,"[","sequence"),recursive = F),nchar),mean)
    )

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
  save(file = 'mergers.rda', mergers, mergers.overlap, mergers.no.overlap, amplicon.info, dadaFs, dadaRs, args)

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

save(out,dadaFs,dadaRs,mergers,seqtab,seqtab.nochim, file = args$dada2_rdata_output)
