library(BSgenome)
library(Biostrings)
library(stringr)
library(dplyr)

args = commandArgs(trailingOnly=T)
numargs=length(args)
ampliconFILE=args[numargs - 2]
genome=args[numargs - 1]
output=args[numargs]

amplicon_info <- read.table(ampliconFILE, header = TRUE)
ref_sequences <- Biostrings::readDNAStringSet(genome)
final_seqs <- NULL

for (idx in 1:nrow(amplicon_info)) {
  info <- amplicon_info[idx, ]
  split_info <- strsplit(info[["amplicon"]], "-")
  start <- info[["ampInsert_start"]]
  end <- info[["ampInsert_end"]]
  pool <- unlist(split_info)[4]
  chr <- unlist(split_info)[1]

  # 'rs' is the reference amplicon sequence
  # 's' is the sequence to make rs
  s <- ref_sequences[str_detect(names(ref_sequences), chr), ]
  if (length(s) == 0) {
    print(paste("skipping", chr))
    next
  }

  rs <- Biostrings::subseq(s, start = start + 1, end = end - 1)

  names(rs) <- paste(c(chr,
                       info[["amplicon_start"]],
                       info[["amplicon_end"]],
                       pool), collapse = "-")

  final_seqs <- c(final_seqs, as.character(rs[1]))
}

set <- DNAStringSet(final_seqs)
Biostrings::writeXStringSet(set, output)