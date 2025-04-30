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
  start <- info[["insert_start"]]
  end <- info[["insert_end"]]
  chr <- info[["chrom"]]

  # 'rs' is the reference amplicon sequence
  # 's' is the sequence to make rs
  s <- ref_sequences[str_detect(names(ref_sequences), chr), ]

  # Error if chromosome not in genome file
  if (length(s) == 0) {
    stop(paste0(
      "ERROR: Chromosome '", chr, "' not found in genome file. Please check if the chromosome names match."
    ))
  }

  rs <- Biostrings::subseq(s, start = start + 2, end = end - 1) # TODO: Add in option to say how much you want to trim left and right by 

  names(rs) <- names(rs) <- info[["target_id"]]

  as.character(rs[1])
}

set <- DNAStringSet(final_seqs)
Biostrings::writeXStringSet(set, args$output)