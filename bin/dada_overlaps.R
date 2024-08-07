library(tidyverse)
library(dada2)
library(argparse)

parser <- ArgumentParser(prog = "DADA2 Workflow", description = "DADA2 workflow script")
parser$add_argument("--trimmed-path",
  type = "character",
  help = "homopolymer threshold to begin masking", nargs = "+", required = TRUE
)
parser$add_argument("--ampliconFILE", type = "character", required = TRUE)
parser$add_argument("--pool", type = "character", default = "pseudo")
parser$add_argument("--band-size", type = "integer", default = 16)
parser$add_argument("--omega-a", type = "double", default = 1e-120)
parser$add_argument("--concat-non-overlaps", action = "store_true")
parser$add_argument("--use-quals", type = "character", default = "false")
parser$add_argument("--maxEE", type = "integer", default = 3)
parser$add_argument("--self-consist", action = "store_true")
parser$add_argument("--omega-c", type = "double", default = 1e-40)
parser$add_argument("--cores", type = "integer", default = 1)


args <- parser$parse_args()
print(args)

fnFs <- sort(list.files(path = args$trimmed_path, pattern = "_R1.fastq.gz", recursive = T, full.names = TRUE))
fnRs <- sort(list.files(path = args$trimmed_path, pattern = "_R2.fastq.gz", recursive = T, full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
print(sample.names)

filtFs <- paste0(args$trimmed_path, "/filtered/", sample.names, "_F_filt.fastq.gz")
filtRs <- paste0(args$trimmed_path, "/filtered/", sample.names, "_R_filt.fastq.gz")


names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  maxN = 0, maxEE = c(args$maxEE, args$maxEE), truncQ = c(5, 5), rm.phix = TRUE,
  compress = TRUE, multithread = args$cores,
  trimRight = c(0, 0), trimLeft = 1, minLen = 75, matchIDs = TRUE
)

filtered_Fs <- filtFs[out[, 2] > 0]
filtered_Rs <- filtRs[out[, 2] > 0]

errF <- learnErrors(filtered_Fs, multithread = args$cores, MAX_CONSIST = 10, randomize = TRUE)
errR <- learnErrors(filtered_Rs, multithread = args$cores, MAX_CONSIST = 10, randomize = TRUE)


pool <- switch(args$pool,
  "true" = TRUE,
  "false" = FALSE,
  "pseudo" = "pseudo"
)

dadaFs <- dada(filtered_Fs, err = errF, selfConsist = args$self_consist, multithread = args$cores, verbose = FALSE, pool = pool, BAND_SIZE = args$band_size, OMEGA_A = args$omega_a)
dadaRs <- dada(filtered_Rs, err = errR, selfConsist = args$self_consist, multithread = args$cores, verbose = FALSE, pool = pool, BAND_SIZE = args$band_size, OMEGA_A = args$omega_a)

if (args$concat_non_overlaps) {
  amplicon.info <- data.frame(
    names = sapply(strsplit(names(dadaFs), "_trimmed"), "[", 1)
  )

  amplicon.info <- amplicon.info %>%
    mutate(names = sapply(str_split(names, "_S(\\d+)"), head, 1)) %>%
    mutate(amplicon = unlist(lapply(str_split(names, "_"), function(x) {
      paste(x[1:3], collapse = "_")
    }))) %>%
    inner_join(
      read.table(args$ampliconFILE, header = T) %>%
        select(amplicon, ampInsert_length),
      by = c("amplicon")
    ) %>%
    select(amplicon, ampInsert_length) %>%
    mutate(
      sum.mean.length.reads = sapply(sapply(unlist(lapply(dadaFs, "[", "sequence"), recursive = F), nchar), mean) +
        sapply(sapply(unlist(lapply(dadaRs, "[", "sequence"), recursive = F), nchar), mean)
    )

  rownames(amplicon.info) <- names(dadaFs)

  mergers.overlap <- mergePairs(dadaFs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length + 10) < sum.mean.length.reads))],
    filtered_Fs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length + 10) < sum.mean.length.reads))],
    dadaRs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length + 10) < sum.mean.length.reads))],
    filtered_Rs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length + 10) < sum.mean.length.reads))],
    verbose = TRUE,
    justConcatenate = FALSE,
    trimOverhang = TRUE,
    minOverlap = 10,
    maxMismatch = 1
  )

  mergers.no.overlap <- mergePairs(dadaFs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length + 10) >= sum.mean.length.reads))],
    filtered_Fs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length + 10) >= sum.mean.length.reads))],
    dadaRs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length + 10) >= sum.mean.length.reads))],
    filtered_Rs[names(dadaFs) %in% rownames(amplicon.info %>% filter((ampInsert_length + 10) >= sum.mean.length.reads))],
    verbose = TRUE,
    justConcatenate = TRUE
  )


  # the dada2 documentation says that a list of data.frames are returned by mergePairs if a list is
  # provided to it, but it also says it can return a data.frame or a list of data.frames. It seems
  # to return a single data.frame if there was only one demultiplex amplicon that did not merge, so
  # convert the output to a list of data.frames so that it can be concatenated with the other merged sequences
  if (class(mergers.no.overlap) == "data.frame") {
    mergers.no.overlap <- list(mergers.no.overlap)
  }

  mergers <- c(mergers.overlap, mergers.no.overlap)[names(dadaFs)]
  save(file = "mergers.rda", mergers, mergers.overlap, mergers.no.overlap, amplicon.info, dadaFs, dadaRs, args)
} else {
  mergers <- mergePairs(dadaFs,
    filtered_Fs,
    dadaRs,
    filtered_Rs,
    verbose = TRUE,
    justConcatenate = FALSE,
    trimOverhang = TRUE,
    minOverlap = 10,
    maxMismatch = 1
  )
}

rm(dadaFs, dadaRs)

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

### Create clusters output file with quality control

seqtab.nochim.df <- as.data.frame(seqtab.nochim)
seqtab.nochim.df$sample <- rownames(seqtab.nochim)

rm(seqtab.nochim)

seqtab.nochim.df[seqtab.nochim.df == 0] <- NA

amplicon.table <- read.table(args$ampliconFILE, header = TRUE, sep = "\t")

# find amplicons (use to select)
pat <- paste(amplicon.table$amplicon, collapse = "|") # find amplicons specified in amplicon table

# create regex to extract sampleID (to be used with remove)
pat.sampleID <- paste(sprintf("^%s_", amplicon.table$amplicon), collapse = "|") # find amplicons specified in amplicon table (with _)
pat.sampleID <- paste(c(pat.sampleID, "_trimmed_merged.fastq.gz$"), collapse = "|") # remove _trimmed from end

seqtab.nochim.df <- seqtab.nochim.df %>%
  pivot_longer(cols = seq(1, ncol(seqtab.nochim.df) - 1), names_to = "asv", values_to = "reads", values_drop_na = TRUE) %>%
  mutate(locus = str_extract(sample, pat)) %>%
  mutate(sampleID = str_remove_all(sample, pat.sampleID)) %>%
  select(sampleID, locus, asv, reads)

allele_sequences <- seqtab.nochim.df %>%
  dplyr::select(locus, asv) %>%
  dplyr::distinct() %>%
  dplyr::group_by(locus) %>%
  dplyr::mutate(allele = seq(1, dplyr::n())) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(allele = paste0(locus, ".", allele)) %>%
  dplyr::rename(sequence = asv)

allele.data <- seqtab.nochim.df %>%
  mutate(sampleID = str_remove_all(sampleID, pat = "_trimmed")) %>%
  left_join(allele_sequences %>% select(-locus), by = c("asv" = "sequence")) %>%
  group_by(sampleID, locus) %>%
  mutate(norm.reads.locus = reads / sum(reads)) %>%
  mutate(n.alleles = n()) %>%
  ungroup()

write.table(allele.data, file = "dada2.clusters.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
