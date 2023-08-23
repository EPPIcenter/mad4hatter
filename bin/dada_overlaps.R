library(tidyverse)
library(dada2)
library(argparse)

parser <- ArgumentParser(prog="DADA2 Workflow", description='DADA2 workflow script')
parser$add_argument('--trimmed-path', type="character",
                   help='homopolymer threshold to begin masking', nargs='+', required = TRUE)
parser$add_argument('--ampliconFILE', type="character", required = TRUE)
parser$add_argument('--pool', type="character", default="false")
parser$add_argument('--band-size', type='integer', default=16)
parser$add_argument('--omega-a', type='double', default=1e-120)
parser$add_argument('--concat-non-overlaps', action='store_true')
parser$add_argument('--use-quals', type="character", default="false")
parser$add_argument('--maxEE', type="integer", default=2)
parser$add_argument('--self-consist', action='store_true')
parser$add_argument('--omega-c', type='double', default=1e-40)
parser$add_argument('--cores', type="integer", default=1)


args <- parser$parse_args()
print(args)

# setwd("~/Documents/GitHub/mad4hatter/work/70/97ca171c7057f885c5372d6c2fdc6e")
# args=list()
# args$trimmed_path = list.files(pattern="trimmed_demuxed*",full.names = T)
# args$dada2_rdada_output = "dada2.seqtab.RDS"
# args$pool = "pseudo"
# args$band_size = -1
# args$omega_a = 1e-120
# args$concat_non_overlaps = T
# args$use_quals="false"
# args$maxEE = 2
# args$ampliconFILE = "v4_amplicon_info.tsv"
# args$omega_c=1e-40
# args$self_consist=T


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
              maxN=0, maxEE=c(args$maxEE, args$maxEE), truncQ=c(5,5), rm.phix=TRUE,
              compress=TRUE, multithread=TRUE,
              trimRight = c(0,0),trimLeft = 1, minLen=75,matchIDs=TRUE) 
head(out)


errF <- learnErrors(filtFs[out[,2]>0], multithread=TRUE,MAX_CONSIST=10,randomize=TRUE)
errR <- learnErrors(filtRs[out[,2]>0], multithread=TRUE,MAX_CONSIST=10,randomize=TRUE)

derepFs <- derepFastq(filtFs[out[,2]>0], verbose = TRUE)
derepRs <- derepFastq(filtRs[out[,2]>0], verbose = TRUE)


pool=switch(
  args$pool,
  "true" = TRUE,
  "false" = FALSE,
  "pseudo" = "pseudo"
)

dadaFs <- dada(derepFs, err=errF, selfConsist=args$self_consist, multithread=args$cores, verbose=FALSE, pool=pool, BAND_SIZE=args$band_size, OMEGA_A=args$omega_a)
dadaRs <- dada(derepRs, err=errR, selfConsist=args$self_consist, multithread=args$cores, verbose=FALSE, pool=pool, BAND_SIZE=args$band_size, OMEGA_A=args$omega_a)

if(args$concat_non_overlaps){
  

  amplicon.info = data.frame(
    names = sapply(strsplit(names(dadaFs),"_trimmed"),"[",1)
  )

  amplicon.info <- amplicon.info %>%
    mutate(names=sapply(str_split(names,'_S(\\d+)'),head,1)) %>% 
    mutate(amplicon=unlist(lapply(str_split(names,'_'), function(x) { paste(x[1:3], collapse = "_") }))) %>%
    inner_join(read.table(args$ampliconFILE,header=T) %>% 
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
  
  mergers.overlap = mergePairs(dadaFs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length+10)<sum.mean.length.reads))],
                               derepFs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length+10)<sum.mean.length.reads))], 
                               dadaRs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length+10)<sum.mean.length.reads))], 
                               derepRs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length+10)<sum.mean.length.reads))], 
                               verbose=TRUE, 
                               justConcatenate=FALSE, 
                               trimOverhang = TRUE,
                               minOverlap=10,
                               maxMismatch=1)
  
  mergers.no.overlap = mergePairs(dadaFs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length+10)>=sum.mean.length.reads))],
                                  derepFs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length+10)>=sum.mean.length.reads))], 
                                  dadaRs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length+10)>=sum.mean.length.reads))], 
                                  derepRs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length+10)>=sum.mean.length.reads))], 
                                  verbose=TRUE, 
                                  justConcatenate=TRUE)
  

  # the dada2 documentation says that a list of data.frames are returned by mergePairs if a list is 
  # provided to it, but it also says it can return a data.frame or a list of data.frames. It seems 
  # to return a single data.frame if there was only one demultiplex amplicon that did not merge, so 
  # convert the output to a list of data.frames so that it can be concatenated with the other merged sequences
  if (class(mergers.no.overlap) == "data.frame") {
    mergers.no.overlap = list(mergers.no.overlap)
  }


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

### Create clusters output file with quality control

seqtab.nochim.df = as.data.frame(seqtab.nochim)
seqtab.nochim.df$sample = rownames(seqtab.nochim)
seqtab.nochim.df[seqtab.nochim.df==0]=NA

amplicon.table=read.table(args$ampliconFILE, header = TRUE, sep = "\t")

# find amplicons (use to select)
pat=paste(amplicon.table$amplicon, collapse="|") # find amplicons specified in amplicon table

# create regex to extract sampleID (to be used with remove)
pat.sampleID=paste(sprintf("^%s_", amplicon.table$amplicon), collapse="|") # find amplicons specified in amplicon table (with _)
pat.sampleID=paste(c(pat.sampleID, "_trimmed_merged.fastq.gz$"),collapse="|")  # remove _trimmed from end

seqtab.nochim.df = seqtab.nochim.df %>%
  pivot_longer(cols = seq(1,ncol(seqtab.nochim)),names_to = "asv",values_to = "reads",values_drop_na=TRUE) %>%
  mutate(locus = str_extract(sample, pat)) %>%
  mutate(sampleID = str_remove_all(sample, pat.sampleID)) %>%
  select(sampleID,locus,asv,reads)

temp = seqtab.nochim.df %>% select(locus,asv) %>% distinct()
loci =unique(temp$locus)
k=1
allele.sequences = data.frame(locus = seq(1,nrow(temp)),allele = seq(1,nrow(temp)),sequence = seq(1,nrow(temp)))
for(i in seq(1,length(loci))){
  temp2 = temp %>% filter(locus==loci[i])
  for(j in seq(1,nrow(temp2))){
    allele.sequences$locus[k+j-1] = loci[i]
    allele.sequences$allele[k+j-1] = paste0(loci[i],".",j)
    allele.sequences$sequence[k+j-1] = temp2$asv[j]
  }
  k=k+nrow(temp2)
}

allele.data = seqtab.nochim.df %>%
  mutate(sampleID = str_remove_all(sampleID, pat = "_trimmed")) %>%
  left_join(allele.sequences %>% select(-locus),by=c("asv"="sequence")) %>%
  group_by(sampleID,locus) %>%
  mutate(norm.reads.locus = reads/sum(reads))%>%
  mutate(n.alleles = n()) %>%
  ungroup()

write.table(allele.data,file="dada2.clusters.txt",quote=F,sep="\t",col.names=T,row.names=F)
