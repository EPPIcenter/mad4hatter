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
parser$add_argument('--use-quals', type="character", default="false")
parser$add_argument('--maxEE', type="integer", default=2)
parser$add_argument('--self-consist', action='store_true')
parser$add_argument('--omega-c', type='double', default=1e-40)

args <- parser$parse_args()
print(args)
## For debugging
# args <- list()
# setwd("/home/bpalmer/Documents/work/7d/875b8c404f348200b6eb70f92e3d80")
# args$trimmed_path = list.files(pattern="trimmed_demuxed")
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
# args$args$dada2_rdata_output="seqtab.nochim.RDS"


# args <- parser$parse_args()
# print(args)

fnFs <- sort(list.files(path=args$trimmed_path, pattern="_R1.fastq.gz", recursive=T, full.names = TRUE))
fnRs <- sort(list.files(path=args$trimmed_path, pattern="_R2.fastq.gz", recursive=T, full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
print(sample.names)

filtFs <- paste0(args$trimmed_path,"/filtered/",sample.names,"_F_filt.fastq.gz")
filtRs <- paste0(args$trimmed_path,"/filtered/",sample.names,"_R_filt.fastq.gz")

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
              maxN=0, maxEE=c(args$maxEE, args$maxEE), truncQ=c(5,5), rm.phix=TRUE,
              compress=TRUE, multithread=TRUE,
              trimRight = c(0,0),trimLeft = 1, minLen=75,matchIDs=TRUE)
head(out)
sample.names.1<-sapply(strsplit(rownames(out[out[,2]>0,]), "_R1"), `[`, 1)

for (sampleID in sample.names.1) {
  cat("Merging:", sampleID, "\n")

  fR <- filtFs[grepl(sampleID, filtFs)]
  rR <- filtRs[grepl(sampleID, filtRs)]

  if (file.exists(fR) && file.exists(rR)) {
    # TODO: add comment here describing literals
    system2("flash", args=sprintf("-t %d -f %d -r %d -s %d --quiet --compress %s %s", 4, 150, 252, 17, fR, rR))
    if (file.size("out.hist") > 0) { # this will be zero if nothing merged
      file.rename("out.extendedFrags.fastq.gz", sprintf("%s_merged.fastq.gz", sampleID))
    }

    if (file.exists('out.notCombined_1.fastq.gz') && file.exists('out.notCombined_1.fastq.gz')) {
      file.stats.1=sprintf("%s-1.txt", sampleID)
      file.stats.2=sprintf("%s-2.txt", sampleID)
      system2("gzip", "-l out.notCombined_1.fastq.gz", stdout = file.stats.1)
      system2("gzip", "-l out.notCombined_2.fastq.gz", stdout = file.stats.2)

      file.stats.table.1=read.table(file.stats.1,header = T)
      file.stats.table.2=read.table(file.stats.2,header = T)
      if (file.stats.table.1$uncompressed > 0 && file.stats.table.2$uncompressed > 0) {
        cat("Unmerged reads detected:", sampleID, "\n")

        file.rename("out.notCombined_1.fastq.gz", sprintf("%s_unmerged1.fastq.gz", sampleID))
        file.rename("out.notCombined_2.fastq.gz", sprintf("%s_unmerged2.fastq.gz", sampleID))
      } else {
        file.remove(c("out.notCombined_1.fastq.gz", "out.notCombined_2.fastq.gz"))
      }
    }
  }
}

merged <- sort(list.files(pattern = "merged.fastq.gz"))
unmerged1 = sort(list.files(pattern = "unmerged1.fastq.gz"))
unmerged2 = sort(list.files(pattern = "unmerged2.fastq.gz"))

all.files=c(merged,unmerged1,unmerged2)

err <- learnErrors(all.files, multithread=TRUE,MAX_CONSIST=10,randomize=TRUE,verbose = 0)

pool=switch(
  args$pool,
  "true" = TRUE,
  "false" = FALSE,
  "pseudo" = FALSE
)

if (!dir.exists("dereps")) {
  dir.create("dereps")
}

clusters<-NULL
derep.fs<-gsub(".fastq.gz", "_derep.rds", all.files)
names(derep.fs)<-all.files

for (sampleID in all.files) {
  tryCatch({
    zz <- file("denoising.step1.txt", open="a")
    cat("Denoising (Step 1):", sampleID, "\n")

    derep <- derepFastq(sampleID, verbose = TRUE)
    clusters[[sampleID]] <- dada(
      derep=derep,
      err=err,
      selfConsist=F,
      multithread=TRUE, verbose=TRUE,
      pool=pool,
      BAND_SIZE=args$band_size,
      OMEGA_A=args$omega_a,
      OMEGA_C=args$omega_c
    )
    saveRDS(derep, file = file.path("dereps", derep.fs[[sampleID]]))
    rm(derep)
  }, error = function(e) {
    message(sprintf("Caught an exception while processing sample %s: %s", sampleID, e$message))
    writeLines(sprintf("%s\t%s\n", sampleID, e$message), con = zz)
  })

  close(zz)
}

## for pseudo pooling - priors are set from the previous run
if ("pseudo" == args$pool) {

  zz <- file("denoising.step2.txt", open="a")
  cat("Denoising with priors...\n")

  seqtab <- makeSequenceTable(clusters)
  clusters <- NULL
  priors <- colnames(seqtab[, colSums(seqtab > 0) >= 2])
  # use the dereplicated sequences (maybe both?)

  for(sampleID in all.files) {
    cat("Denoising (Step 2):", sampleID, "\n")
    tryCatch({
      derep <- readRDS(file.path("dereps", derep.fs[[sampleID]]))
      clusters[[sampleID]] <- dada(
        derep=derep,
        err=err,
        selfConsist=F,
        multithread=TRUE, verbose=TRUE,
        pool=pool,
        BAND_SIZE=args$band_size,
        OMEGA_A=args$omega_a,
        OMEGA_C=args$omega_c
      )

      rm(derep)
    }, error = function(e) {
      message(sprintf("Caught an exception while processing sample %s: %s", sampleID, e$message))
      writeLines(sprintf("%s\t%s\n", sampleID, e$message), con = zz)
    })
  }

  close(zz)
}

unmerged.dereps.1=sort(derep.fs[names(derep.fs) %in% c(unmerged1)])
unmerged.dereps.2=sort(derep.fs[names(derep.fs) %in% c(unmerged2)])
unmerged.merged=list()

for (ii in seq_along(unmerged.dereps.1)) {
  rds.1=file.path("dereps", unmerged.dereps.1[[ii]])
  rds.2=file.path("dereps", unmerged.dereps.2[[ii]])
  if (file.exists(rds.1) && file.exists(rds.2)) {
    derepF=readRDS(rds.1)
    derepR=readRDS(rds.2)
    dadaF=clusters[[names(unmerged.dereps.1[ii])]]
    dadaR=clusters[[names(unmerged.dereps.2[ii])]]
    merged.id=str_replace(names(unmerged.dereps.1[ii]), "unmerged1", "merged")
    unmerged.merged[[merged.id]]=mergePairs(dadaF,derepF,dadaR,derepR,justConcatenate = TRUE)
  } else {
    # print info message here...
  }
}

seqtab.unmerged <- makeSequenceTable(unmerged.merged)
seqtab.merged <- makeSequenceTable(clusters[names(clusters) %in% merged])
seqtab=mergeSequenceTables(seqtab.unmerged,seqtab.merged, repeats = "sum")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

save(clusters, err, seqtab.nochim, file = "DADA2.RData")
saveRDS(seqtab.nochim, file = args$dada2_rdata_output)
saveRDS(clusters, file = "clusters.RDS")
