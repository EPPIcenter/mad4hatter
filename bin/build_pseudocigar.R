#' @title Create the pseudoCIGAR string for aligned ASVs
#' 
#' @description This script produces a pseudoCIGAR string using masked or unmasked ASVs based on the input alignments.
#' 
#' @keywords pseudocigar

library(argparse)

parser <- ArgumentParser(description='Create the pseudoCIGAR string using masked or unmasked ASVs')
parser$add_argument('--alignments', type="character", required=TRUE, help = "File containing aligned ASVs")

args <- parser$parse_args()
print(args)

library(stringr)
library(dplyr)
library(magrittr)

#' Compute Insertion Group
#' 
#' @description Identify the group of insertion for given reference characters.
#' @param ref_chars A vector of characters representing the reference.
#' @return A vector of integers representing insertion groups.
compute_insertion_group <- function(ref_chars) {
  insertion_start <- ifelse(lag(ref_chars, default = "") != "-" & ref_chars == "-", 1, 0)
  return(cumsum(insertion_start))
}

#' Compute Deletion Group
#' 
#' @description Identify the group of deletion for given query characters.
#' @param query_chars A vector of characters representing the query.
#' @return A vector of integers representing deletion groups.
compute_deletion_group <- function(query_chars) {
  deletion_start <- ifelse(lag(query_chars, default = "") != "-" & query_chars == "-", 1, 0)
  return(cumsum(deletion_start))
}

#' Compute Insertion CIGAR String
#' 
#' @description Generate the CIGAR string representation for insertions.
#' @param position Position of insertion.
#' @param query_char Query characters.
#' @param insertion_group Group identifier for insertion.
#' @return A string representing the CIGAR format for insertion.
compute_insertion_cigar <- function(position, query_char, insertion_group) {
  # Get the starting position of the insertion
  start_position <- first(position[insertion_group == unique(insertion_group)[insertion_group]])
  
  # Get the inserted sequence
  inserted_seq <- paste0(query_char[insertion_group == unique(insertion_group)[insertion_group]], collapse = "")
  
  return(paste0(start_position, "I=", inserted_seq))
}

#' Compute Deletion CIGAR String
#' 
#' @description Generate the CIGAR string representation for deletions.
#' @param position Position of deletion.
#' @param ref_char Reference characters.
#' @param deletion_group Group identifier for deletion.
#' @return A string representing the CIGAR format for deletion.
compute_deletion_cigar <- function(position, ref_char, deletion_group) {
  # Get the starting position of the deletion
  start_position <- first(position[deletion_group == unique(deletion_group)[deletion_group]])
  
  # Get the deleted sequence
  deleted_seq <- paste0(ref_char[deletion_group == unique(deletion_group)[deletion_group]], collapse = "")
  
  return(paste0(start_position, "D=", deleted_seq))
}

#' Build PseudoCIGAR String
#' 
#' @description Create a pseudoCIGAR string based on the reference and query sequences.
#' @param reference Reference sequence.
#' @param query Query sequence.
#' @return A string representing the pseudoCIGAR format.
build_pseudoCIGAR_string <- function(reference, query) {

  # Transform strings into vectors of characters
  ref_chars <- str_split(reference, pattern = "")[[1]]
  query_chars <- str_split(query, pattern = "")[[1]]
  
  # Prepare a dataframe
  df <- tibble(position = seq_along(ref_chars),
               ref_char = ref_chars,
               query_char = query_chars)
  
  # Get positions that are in N-, NNNN, or N-N
  n_string <- str_locate_all(reference, pattern = "N+\\-+N*|N+")[[1]]
  n_vector <- unlist(mapply(seq, from = n_string[, "start"], to = n_string[, "end"]))

  # Indicate these positions are included in mask when building cigar
  df %<>% 
    mutate(
      masked = position %in% n_vector,
      insertion_group = compute_insertion_group(ref_chars),
      deletion_group = compute_deletion_group(query_chars)
    ) %>%
    mutate(result = case_when(
      masked ~ paste0(n_string[, "start"], "N+", n_string[, "end"] - n_string[, "start"] + 1),
      ref_char != query_char & ref_char != "-" & query_char != "-" ~ compute_deletion_cigar(position, ref_char, deletion_group), 
      ref_char != query_char & ref_char == "-" & query_char != "-" ~ compute_insertion_group(position, query_char, insertion_group),
      ref_char != query_char & ref_char != "-" & query_char == "-" ~ paste0(position, "D=", ref_char),
      TRUE ~ ""
    ))
  
  # Combine the results into a single string
  cigar_str <- str_c(df$result, collapse = "")
  if (str_length(cigar_str) == 0) {
    cigar_str <- "." # if the cigar string is empty, replace with a dot
  }

  return (cigar_str)
}

# Main code to create the pseudocigar string

# Load in alignment data. This can be masked (with 'N's)or unmasked data.
df.aln = read.csv(args$alignments, sep="\t", header=T)

# Create the pseudoCIGAR string fror each alignment entry
pseudo_cigar = df.aln %>%
  dplyr::mutate(pseudo_cigar = mapply(build_pseudoCIGAR_string, refseq, hapseq))

# Write out the data frame
write.table(pseudo_cigar,file="alignments.pseudocigar.txt",quote=F,sep="\t",col.names=T,row.names=F)
