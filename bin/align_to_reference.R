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

# Example of usage
load_library("dplyr")
load_library("Biostrings")

load_library("argparse")
parser <- ArgumentParser(description='Aligns sequences (DADA2 clusters) against a reference of known sequences')
parser$add_argument('--clusters', type="character", required=TRUE, help="RDS Clusters from DADA2. This is the main output from the DADA module.")
parser$add_argument('--refseq-fasta', type="character")
parser$add_argument('--alignment-threshold', type="integer", default = 60)
parser$add_argument('--n-cores', type = 'integer', default = 1, help = "Number of cores to use. Ignored if running parallel flag is unset.")
parser$add_argument('--amplicon-table', type="character", required=TRUE, help = "Amplicon table with primer pools. This is used to organize the sequence table by amplicon.")
parser$add_argument('--log-level', type="character", default = "INFO", help = "Log level. Default is INFO.")

args <- parser$parse_args()

log_level_arg <- match.arg(args$log_level, c("DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
log_threshold(log_level_arg)
log_appender(appender_console)
args_string <- paste(sapply(names(args), function(name) {
  paste(name, ":", args[[name]])
}), collapse = ", ")
log_debug(args_string)

load_library("stringr")
load_library("dplyr")
load_library("dada2")
load_library("foreach")
load_library("parallel")
load_library("BSgenome")
load_library("tidyr")
load_library("doMC")
load_library("tibble")
load_library("ggplot2")
load_library("Biostrings")
load_library("magrittr")

## FOR DEUBGING
# setwd("/home/bpalmer/Documents/GitHub/mad4hatter/work/48/4bb9040a0aa3b2febe32ce38f4e3bd")
# args=list()
# args$clusters="clusters.concatenated.collapsed.txt"
# args$refseq_fasta="v4_reference.fasta"
# args$parallel=TRUE
# args$n_cores=2
# args$amplicon_table="v4_amplicon_table.tsv"

clusters=read.table(args$clusters, header=T)
log_info("Clusters read from: {args$clusters}")

# register number of cores to use
registerDoMC(args$n_cores)
log_info("Registered {args$n_cores} cores for parallel processing")

# read the sequences from the reference genome (already extracted into a fasta file)
ref_sequences <- readDNAStringSet(args$refseq_fasta)
log_info("Reference sequences read from: {args$refseq_fasta}")

# DADA2 should have corrected any substitution errors during sequencing. Therefore, any 
# mismatches are expected to be real. A small penalty is added when there is a mismatching base
# to catch large differences that are likely not real. Gap penalties are also set below to 
# help filter sequences that are truly different from the reference. 
match <- 2
mismatch <- -1
baseOnly <- TRUE
sigma <- nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = baseOnly)
log_info("Substitution matrix initialized: match = {match}, mismatch = {mismatch}, baseOnly = {baseOnly}")

#########  IMPORTANT: NON PF SEQUENCES NEED TO BE INCLUDED IN REFERENCE!

# This object contains the aligned ASV sequences
df_aln <- foreach(seq1 = 1:nrow(clusters), .combine = "bind_rows") %dopar% {
  log_debug("Processing sequence {seq1} of {nrow(clusters)}")
  # the alignment is performed only with the reference sequence for the corresponding locus as we have this information from the demultiplexing step
  refseq.seq1 = ref_sequences[clusters$locus[seq1]] #commenting out next lines as I concatenated all genomes
  aln <- pairwiseAlignment(refseq.seq1, str_remove_all(clusters$asv[seq1],"N"), substitutionMatrix = sigma, gapOpening = -8, gapExtension = -5, scoreOnly = FALSE)
  patt <- c(alignedPattern(aln), alignedSubject(aln))
  ind <- sum(str_count(as.character(patt),"-"))
  data.frame(
    sampleID = clusters$sampleID[seq1],
    asv = clusters$asv[seq1],
    hapseq = as.character(patt)[2],
    refseq = as.character(patt)[1],
    refid = clusters$locus[seq1],
    score = score(aln),
    indels = ind
  )
}


write.table(df_aln,file="alignments.txt",quote=F,sep="\t",col.names=T,row.names=F)
log_info("Alignments written to alignments.txt")
