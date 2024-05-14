library(BSgenome)
library(Biostrings)
library(stringr)
library(dplyr)
library(argparse)
library(foreach)
library(doMC)

parser <- ArgumentParser(description='Create reference sequences using amplicon table')
parser$add_argument('--output', type="character", help='name of fasta to output', required = TRUE)
parser$add_argument('--ampliconFILE', type="character", required = TRUE)
parser$add_argument('--genome', type="character", required = TRUE)
parser$add_argument('--ncores', type="integer", default=1)

args <- parser$parse_args()

amplicon_info <- read.table(args$ampliconFILE, header = TRUE)
ref_sequences <- Biostrings::readDNAStringSet(args$genome)

doMC::registerDoMC(cores = args$ncores)
final_seqs <- foreach (idx = 1:nrow(amplicon_info), .combine = "c") %dopar% {
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
    stop("Missing reference detected!!!")
  }

  # The positioning in the amplicon panel file is 0-based.
  # Thus, if a `start_position` is 1397996, then the actual
  # position in the sequence is 1397997. The start_position is
  # the actual position where the amplicon begins, and starts with the
  # primer sequence. The length of the primer sequences is added to
  # identify the amplicon insert position. The position of the actual
  # ASV is one base more, because it starts after the primer sequence.
  rs <- Biostrings::subseq(s, start = start + 1, end = end - 1)

  names(rs) <- paste(c(chr,
                       info[["amplicon_start"]],
                       info[["amplicon_end"]],
                       pool), collapse = "-")

  as.character(rs[1])
}

set <- DNAStringSet(final_seqs)
Biostrings::writeXStringSet(set, args$output)