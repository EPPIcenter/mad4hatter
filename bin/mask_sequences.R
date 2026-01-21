# Import necessary libraries
library(argparse)
library(Biostrings)
library(dplyr)
library(purrr)
library(doMC)
library(parallel)
library(stringr)
library(magrittr)

# ---------------------------
# DEFINE ARGUMENT PARSER
# ---------------------------

# Description of the script
parser <- ArgumentParser(description = "Post process dada2 inferred sequences")

# Define input arguments
parser$add_argument("--alignments", type = "character", required = TRUE, 
                    help = "RDS Clusters from DADA2. Main output from DADA module.")
parser$add_argument("--masks", type = "character", nargs = "+")
parser$add_argument("--parallel", action = "store_true")
parser$add_argument("--n-cores", type = "integer", default = -1, 
                    help = "Number of cores to use. Ignored if running parallel flag is unset.")

# Parse the arguments
args <- parser$parse_args()
print(args)

# Read alignment data
df_aln <- read.table(args$alignments, header = TRUE)

# ---------------------------
# SET UP PARALLEL PROCESSING
# ---------------------------

# Setup parallel backend if asked
if (args$parallel) {
  n_cores <- ifelse(args$n_cores <= 0, detectCores(), args$n_cores)
  registerDoMC(n_cores)
} else {
  registerDoSEQ()
}

# ---------------------------
# DEFINE FUNCTIONS
# ---------------------------

# Function to find 'N' positions
get_masked_pos <- function(sequence) {
  out <- unlist(gregexpr("N", sequence, ignore.case = TRUE))
  if (length(out) == 1 && out == -1) {
    return (NA)
  } else {
    return (out)
  }
}

# Function to get positions of 'N' from a FASTA file
get_fasta_masked_positions <- function(fasta_path) {
  fasta <- readDNAStringSet(fasta_path)
  positions_df <- data.frame(refid = names(fasta), stringsAsFactors = FALSE)
  positions_df$pos <- lapply(fasta, function(x) get_masked_pos(as.character(x)))
  names(positions_df)[2] <- fasta_path 
  return(positions_df)
}

# Function to mask aligned sequence
mask_aligned_sequence <- function(pos, refseq) {
  if (is.null(pos) || (is.list(pos) && length(pos) == 0)) {
    return (refseq)
  }
  pattern <- gregexpr("[^-]+|-+", refseq)
  split_string <- regmatches(refseq, pattern)[[1]]
  string_homo_tr_ref_seq <- str_remove_all(refseq, "-")
  string_homo_tr_ref_seq <- strsplit(as.character(string_homo_tr_ref_seq), "")[[1]]
  string_homo_tr_ref_seq[unlist(pos)] <- "N"
  string_homo_tr_ref_seq <- paste(string_homo_tr_ref_seq, collapse = "")
  newstring <- character()
  for(chunk in split_string) {
    if(grepl("-", chunk)) {
      newstring <- paste0(newstring, chunk)
    } else {
      newstring <- paste0(newstring, substr(string_homo_tr_ref_seq, 1, nchar(chunk)))
      string_homo_tr_ref_seq <- substr(string_homo_tr_ref_seq, nchar(chunk) + 1,
                                       nchar(string_homo_tr_ref_seq))
    }
  }
  return(newstring)
}

# ---------------------------
# EXECUTE MAIN SCRIPT
# ---------------------------

# Apply the function to get masked positions
fasta_positions <- lapply(args$masks, get_fasta_masked_positions)

# Merge all data frames
merged <- Reduce(function(df1, df2) merge(df1, df2, by = "refid", all = TRUE), fasta_positions)
merged$all_positions <- apply(merged[-which(names(merged) %in% "refid")], 1,
                               function(x) unlist(x[!is.na(unlist(x))]))
merged <- merged[, c("refid", "all_positions")]

# Mask the aligned reference sequences
df_aln_references <- df_aln %>% select(refid, refseq) %>% distinct() %>%
                     left_join(merged, by = "refid")

df_aln_references %<>% 
  mutate(masked_aligned_refseq = mapply(mask_aligned_sequence, all_positions, refseq))

# Mask the hapseq sequences 
df_aln_hapseq <- df_aln %>% select(refid,hapseq) %>% distinct() %>%
  left_join(merged, by = "refid")

df_aln_hapseq %<>% 
  mutate(masked_hapseq = mapply(mask_aligned_sequence, all_positions, hapseq))

# Update df_aln
df_aln %<>% left_join(df_aln_references, by = c("refid", "refseq")) %>%
  left_join(df_aln_hapseq, by = c("refid", "hapseq")) %>%
  select(sampleID, asv, refid, hapseq, masked_aligned_refseq, masked_hapseq, score, indels) %>%
  dplyr::rename(refseq = masked_aligned_refseq)

# Write the output to a file
write.table(df_aln, file = "masked.alignments.txt", quote = FALSE, sep = "\t", 
            col.names = TRUE, row.names = FALSE)
