
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


# Brian Code
v <- vector(length = length(mergers))
for (idx in seq_along(mergers)) {
    v[idx] <- nchar(mergers[[idx]]$sequence) > 0 &&
        isTRUE(str_detect(mergers[[idx]]$sequence, "NNNNNNNNNN"))
}

merged_non_overlappers <- mergers[v]

aln_results <- list()
if (length(merged_non_overlappers) > 0) {

    # get the substitution scores given a paremterized error probability
    subscores <- errorSubstitutionMatrices(
        errorProbability = 0.1) # parameterize this

    submat <- matrix(subscores[, , "0"], 4, 4)
    diag(submat) <- subscores[, , "1"]
    dimnames(submat) <- list(DNA_ALPHABET[1:4], DNA_ALPHABET[1:4])

    for (merged_non_overlap in merged_non_overlappers) {

        best_score <- -Inf
        best_alignment <- NULL

        # get the starting sequence
        starting_sequence_name <- names(merged_non_overlap)
        starting_sequence <- DNAString(merged_non_overlappers[[
            starting_sequence_name]]$sequence)

        # get the amplicon name from concatenated sample-amplicon name
        amplicon_name <- sub(
            "_PARAV3.*", "", names(merged_non_overlap)) # parameterize this

        # get all amplicons by name across samples
        amplicons_across <- str_detect(
            names(mergers), amplicon_name)

        for (amplicon in mergers[amplicons_across]) {

            # make sure were not comparing against ourselves
            if (!identical(names(amplicon), names(merged_non_overlap))) {
                # Fit a pairwise comparison
                compare_sequence <- DNAString(
                    mergers[[names(amplicon)]]$sequence)

                score <- pairwiseAlignment(starting_sequence, compare_sequence,
                    substitutionMatrix = subscores, scoreOnly = TRUE)

                if (score > best_score) {
                    best_alignment <- amplicon
                    best_score <- score
                }
            }
        }

        aln_results <- c(aln_results, list(
            starting_sequence_name = starting_sequence_name,
            best_alignment = best_alignment,
            best_score = best_score)
        )
    }
}

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

save(out,dadaFs,dadaRs,mergers,seqtab,seqtab.nochim,aln_results, file = dada2RData)
