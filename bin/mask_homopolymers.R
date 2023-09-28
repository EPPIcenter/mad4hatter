library(logger)
load_library <- function(library_name) {
  output <- capture.output({
    suppressWarnings({
      library(library_name, character.only = TRUE)
    })
  }, type = "message")
  
  # Separate warnings from messages
  warnings <- warnings()
  
  # Log messages
  if(length(output) > 0) {
    log_info(paste("Message from", library_name, ":", paste(output, collapse = "; ")))
  }
  
  # Log warnings
  if(length(warnings) > 0) {
    log_warn(paste("Warning from", library_name, ":", paste(warnings, collapse = "; ")))
  }
}

load_library("argparse")

parser <- ArgumentParser(description='Mask homopolymers in a fasta')
parser$add_argument('--refseq-fasta', type="character", required=TRUE)
parser$add_argument('--homopolymer_threshold', type="integer", default = 5)
parser$add_argument('--fout', type="character", default="refseq.homopolymer.fasta.mask")
parser$add_argument('--n-cores', type="integer", default=4)
parser$add_argument('--log-level', type="character", default = "INFO", help = "Log level. Default is INFO.")

args <- parser$parse_args()
args_string <- paste(sapply(names(args), function(name) {
  paste(name, ":", args[[name]])
}), collapse = ", ")

log_debug(paste("Arguments parsed successfully:", args_string))

load_library("BSgenome")
load_library("Biostrings")
load_library("BiocParallel")
load_library("dplyr")

# Define a function to process each sequence
homomask_sequence <- function(sequence) {
  # somewhat convoluted way to get the Rle from the sequence in the DNAStringSet object that will be passed to the function
  ref_seq_rle <- Rle(as.vector(unlist(strsplit(as.character(sequence), ""))))
  # replace the homopolymers with Ns
  ref_seq_rle@values[ref_seq_rle@lengths > args$homopolymer_threshold] <- "N"
  # add 2 N's to the run to mask 1 base upstream and 1 base downstream
  ref_seq_rle@lengths[ref_seq_rle@lengths > args$homopolymer_threshold] <-  as.integer(ref_seq_rle@lengths[ref_seq_rle@lengths > args$homopolymer_threshold] + 2)
  # remove one N from the run to mask 1 base upstream and 1 base downstream
  ref_seq_rle@lengths[lag(ref_seq_rle@lengths > args$homopolymer_threshold,default=FALSE)] <-  as.integer(ref_seq_rle@lengths[lag(ref_seq_rle@lengths > args$homopolymer_threshold,default=FALSE)] -1)
  ref_seq_rle@lengths[lead(ref_seq_rle@lengths > args$homopolymer_threshold,default=FALSE)] <-  as.integer(ref_seq_rle@lengths[lead(ref_seq_rle@lengths > args$homopolymer_threshold,default=FALSE)] -1)
  # if a base was between homopolymers we have a -1 now, so we need to account by taking an N out of one of the homopolymers, and I'll do from the left
  ref_seq_rle@lengths[lead(ref_seq_rle@lengths == -1 ,default=FALSE)] <-  as.integer(ref_seq_rle@lengths[lead(ref_seq_rle@lengths == -1 ,default=FALSE)]  -1)
  # if homopolymer was at the beginning or end we added an N, so let's remove it
  ref_seq_rle@lengths[c(1, length(ref_seq_rle@lengths))] <- ref_seq_rle@lengths[c(1, length(ref_seq_rle@lengths))] - (ref_seq_rle@values[c(1, length(ref_seq_rle@lengths))] == "N")
  # remove empty runs, make -1's, 0
  ref_seq_rle <- Rle(ref_seq_rle@values, ifelse(ref_seq_rle@lengths<0, 0, ref_seq_rle@lengths))
  # paste the Rle back together into a sequence
  homomasked_seq <- DNAString(paste(rep(runValue(ref_seq_rle), times = runLength(ref_seq_rle)), collapse = ""))
  # return the sequence
  return(homomasked_seq)
}

register(MulticoreParam(workers = args$n_cores)) # adjust 'workers' to the number of cores you want to use
inDNAStringSet <- readDNAStringSet(args$refseq_fasta)
homomask_ref_sequences <- DNAStringSet(bplapply(inDNAStringSet, homomask_sequence)) # mask homopolymers
writeXStringSet(homomask_ref_sequences, file=args$fout) # write out