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

load_library("BSgenome")
load_library("Biostrings")
load_library("stringr")
load_library("dplyr")
load_library("argparse")
load_library("foreach")
load_library("doMC")

parser <- ArgumentParser(description='Create reference sequences using amplicon table')
parser$add_argument('--output', type="character", help='name of fasta to output', required = TRUE)
parser$add_argument('--ampliconFILE', type="character", required = TRUE)
parser$add_argument('--genome', type="character", required = TRUE)
parser$add_argument('--ncores', type="integer", default=1)
parser$add_argument('--log-level', type="character", default = "INFO", help = "Log level. Default is INFO.")

args <- parser$parse_args()
# Set up logging
log_threshold(args$log_level)
log_appender(appender_console)
args_string <- paste(sapply(names(args), function(name) {
  paste(name, ":", args[[name]])
}), collapse = ", ")

log_debug(paste("Arguments parsed successfully:", args_string))

amplicon_info <- read.table(args$ampliconFILE, header = TRUE)
ref_sequences <- Biostrings::readDNAStringSet(args$genome)
log_info(paste("Reference sequences read from:", args$genome))

doMC::registerDoMC(cores = args$ncores)
log_info(paste("Using", args$ncores, "cores for parallel processing"))

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
    print(paste("skipping", chr))
    return(NULL)
  }

  rs <- Biostrings::subseq(s, start = start + 1, end = end - 1)

  names(rs) <- paste(c(chr,
                       info[["amplicon_start"]],
                       info[["amplicon_end"]],
                       pool), collapse = "-")

  as.character(rs[1])
}

set <- DNAStringSet(final_seqs)
Biostrings::writeXStringSet(set, args$output)
log_info(paste("Reference sequences written to:", args$output))