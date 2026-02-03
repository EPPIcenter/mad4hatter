library(BSgenome)
library(Biostrings)
library(stringr)
library(dplyr)
library(argparse)
library(foreach)
library(doMC)

parser <- ArgumentParser(description = "Create reference sequences using amplicon table")
parser$add_argument("--output", type = "character", help = "name of fasta to output", required = TRUE)
parser$add_argument("--ampliconFILE", type = "character", required = TRUE)
parser$add_argument("--genome", type = "character", required = TRUE)
parser$add_argument("--ncores", type = "integer", default = 1)

args <- parser$parse_args()

amplicon_info <- read.table(args$ampliconFILE, header = TRUE)
ref_sequences <- Biostrings::readDNAStringSet(args$genome)

# Check that all chromosomes exist in the genome file before processing
unique_chromosomes <- unique(amplicon_info[["chrom"]])
missing_chromosomes <- c()
for (chr in unique_chromosomes) {
  if (length(ref_sequences[str_detect(names(ref_sequences), chr), ]) == 0) {
    missing_chromosomes <- c(missing_chromosomes, chr)
  }
}

if (length(missing_chromosomes) > 0) {
  stop(paste("Error: The following chromosomes were not found in the genome file:",
             paste(missing_chromosomes, collapse = ", ")))
}

doMC::registerDoMC(cores = args$ncores)
final_seqs <- foreach(idx = 1:nrow(amplicon_info), .combine = "c") %dopar% {
  info <- amplicon_info[idx, ]
  start <- info[["insert_start"]]
  end <- info[["insert_end"]]
  chr <- info[["chrom"]]

  # 'rs' is the reference amplicon sequence
  # 's' is the sequence to make rs
  s <- ref_sequences[str_detect(names(ref_sequences), chr), ]
  # This check should never be true now, but kept as a safety measure
  if (length(s) == 0) {
    stop(paste("Error: Chromosome", chr, "not found in genome file"))
  }
  # biostrings takes 1-based coordinates. We want to extract with trim of 1 base from each end
  #TODO: allow trim to be input by user
  rs <- Biostrings::subseq(s, start = start + 2, end = end - 1)

  names(rs) <- info[["target_name"]]

  as.character(rs[1])
}

set <- DNAStringSet(final_seqs)
Biostrings::writeXStringSet(set, args$output)
