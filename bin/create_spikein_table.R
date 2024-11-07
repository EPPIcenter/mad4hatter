library(tidyverse)
library(dada2)
library(argparse)

parser <- ArgumentParser(prog="DADA2 Workflow", description='DADA2 workflow script')
parser$add_argument('--trimmed-path', type="character",
                   help='spikein fastqs to include in spike table', nargs='+', required = TRUE)
parser$add_argument('--maxEE', type="integer", default=2)

args <- parser$parse_args()
print(args)

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
derep <- derepFastq(merged, verbose = TRUE)

seqtab <- makeSequenceTable(derep)
seqtab_df <- as.data.frame(seqtab)
seqtab_df$SampleID <- rownames(seqtab_df)

# find spikein ID (use to select)
pat <- 'SDSI'

# create regex to extract sampleID (to be used with remove)
pat.sampleID <- paste(sprintf("^%s_", pat), collapse = "|") # find amplicons specified in amplicon table (with _)
pat.sampleID <- paste(c(pat.sampleID, "_trimmed_merged.fastq.gz$"), collapse = "|") # remove _trimmed from end

seqtab_df <- seqtab_df %>%
  pivot_longer(cols = seq(1, ncol(seqtab_df) - 1), names_to = "Spikein", values_to = "Reads", values_drop_na = TRUE) %>%
  mutate(SampleID = str_remove_all(SampleID, pat.sampleID)) %>%
  mutate(SampleID = str_remove_all(SampleID, "_merged.fastq.gz$")) %>%
  select(SampleID, Spikein, Reads)

# seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# save(clusters, err, seqtab.nochim, file = "DADA2.RData")
# saveRDS(seqtab.nochim, file = args$dada2_rdata_output)
write.table(seqtab_df, file = "spikeins_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)
