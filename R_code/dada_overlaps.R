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
args <- list()
setwd("/home/bpalmer/Documents/work/5e/0c033ae668caa0c89d4d7cbe574985")
args$trimmed_path = file.path("/home/bpalmer/Documents/work/2c/1c3462361b2882484daec490c7859e", c("trimmed_demuxed1", "trimmed_demuxed2"))
args$dada2_rdada_output = "dada2.seqtab.RDS"
args$pool = "pseudo"
args$band_size = -1
args$omega_a = 1e-120
args$concat_non_overlaps = T
args$use_quals="false"
args$maxEE = 2
args$ampliconFILE = "v4_amplicon_info.tsv"
args$omega_c=1e-40
args$self_consist=T


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
  }
}

merged <- sort(list.files(pattern = "merged.fastq.gz"))
err <- learnErrors(merged, multithread=TRUE,MAX_CONSIST=10,randomize=TRUE,verbose = 0)

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
derep.fs<-gsub(".fastq.gz", "_derep.rds", merged)
names(derep.fs)<-merged

for (sampleID in merged) {
  tryCatch({
    cat("Denoising (Step 1):", sampleID, "\n")

    derep <- derepFastq(sampleID, verbose = TRUE)
    clusters[[sampleID]] <- dada(
      derep=derep,
      err=err,
      selfConsist=args$self_consist,
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
  })
}

## for pseudo pooling - priors are set from the previous run
if ("pseudo" == args$pool) {

  cat("Denoising with priors...\n")

  seqtab <- makeSequenceTable(clusters)
  clusters <- NULL
  priors <- colnames(seqtab[, colSums(seqtab > 0) >= 2])
  # use the dereplicated sequences (maybe both?)

  for(sampleID in merged) {
    cat("Denoising (Step 2):", sampleID, "\n")
    tryCatch({
      derep <- readRDS(file.path("dereps", derep.fs[[sampleID]]))
      clusters[[sampleID]] <- dada(
        derep=derep,
        err=err,
        selfConsist=args$self_consist,
        multithread=TRUE, verbose=TRUE,
        pool=pool,
        BAND_SIZE=args$band_size,
        OMEGA_A=args$omega_a,
        OMEGA_C=args$omega_c
      )

      rm(derep)
    }, error = function(e) {
      message(sprintf("Caught an exception while processing sample %s: %s", sampleID, e$message))
    })
  }
}

seqtab <- makeSequenceTable(clusters)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

save(clusters, err, seqtab.nochim, file = "DADA2.RData")
saveRDS(seqtab.nochim, file = args$dada2_rdata_output)
saveRDS(clusters, file = "clusters.RDS")
