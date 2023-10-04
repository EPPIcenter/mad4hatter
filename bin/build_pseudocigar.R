#' @title Create the pseudoCIGAR string for aligned ASVs
#' 
#' @description This script produces a pseudoCIGAR string using masked or unmasked ASVs based on the input alignments.
#' 
#' @keywords pseudocigar

library(stringr)
library(dplyr)
library(magrittr)
library(foreach)
library(doMC)

#' Compute Insertion Group
#'
#' @description Identify the group of insertion for given reference characters.
#' @param ref_chars A vector of characters representing the reference.
#' @param mask_group A vector of integers representing mask groups.
#' @return A vector of integers representing insertion groups.
compute_insertion_group <- function(ref_chars, mask_group) {
  runs <- rle(ref_chars)

  # Identify runs of "-"
  is_insert <- runs$values == "-"

  # Create group numbers for the runs of insertions
  group_nums <- cumsum(is_insert)
  group_nums[!is_insert] <- NA

  # Expand the group numbers based on the lengths of the runs
  result <- rep(group_nums, times = runs$lengths)

  # If mask group has a value, increment the group number for subsequent insertions
  if (!all(is.na(mask_group))) {
    mask_positions <- which(!is.na(mask_group))
    for (pos in mask_positions) {
      result[pos:length(result)] <- result[pos:length(result)] + 1
    }
  }

  return(result)
}

#' Compute Deletion Group
#'
#' @description Identify the group of deletion for given query characters.
#' @param query_chars A vector of characters representing the query.
#' @param mask_group A vector of integers representing mask groups.
#' @return A vector of integers representing deletion groups.
compute_deletion_group <- function(query_chars, mask_group) {
  runs <- rle(query_chars)

  # Identify runs of "-"
  is_deletion <- runs$values == "-"

  # Create group numbers for the runs of deletions
  group_nums <- cumsum(is_deletion)
  group_nums[!is_deletion] <- NA

  # Expand the group numbers based on the lengths of the runs
  result <- rep(group_nums, times = runs$lengths)

  # If mask group has a value, increment the group number for subsequent deletions
  if (!all(is.na(mask_group))) {
    mask_positions <- which(!is.na(mask_group))
    for (pos in mask_positions) {
      result[pos:length(result)] <- result[pos:length(result)] + 1
    }
  }

  return(result)
}

#' Compute Mask Group
#'
#' This function identifies positions in a reference sequence that should be treated as masked.
#' A sequence is considered "masked" if it has consecutive 'N's possibly interspersed
#' with '-'. For example, the sequences "N-", "N-N", "NN-", and "N--N" are all considered
#' masked sequences. The output provides a vector of integers where each masked region is
#' assigned a unique group identifier.
#'
#' @param ref_chars A vector of characters representing the reference sequence.
#'
#' @return A vector of integers representing mask groups. Positions not part of a mask are assigned NA.
#'
#' @examples
#' compute_mask_group(c("A", "N", "N", "-", "G"))
#' # Returns: NA 2 2 2 NA
#'
#' compute_mask_group(c("A", "N", "N", "-", "G", "N", "-", "N", "T"))
#' # Returns: NA 2 2 2 NA 5 5 5 NA
#'
#' compute_mask_group(c("A", "G", "T"))
#' # Returns: NA NA NA
#'
compute_mask_group <- function(ref_chars) {
  # Convert the vector of characters to a single string for regex operations
  ref_string <- paste0(ref_chars, collapse = "")

  # Locate sequences based on the pattern
  mask_regions <- str_locate_all(ref_string, pattern = "N+(-*N*)*")[[1]]

  # Initialize the mask_groups vector with NAs
  mask_groups <- rep(NA, length(ref_chars))

  # For each mask start, fill in all positions up to the mask end
  for (i in seq_along(mask_regions[, "start"])) {
    mask_groups[mask_regions[i, "start"]:mask_regions[i, "end"]] <- mask_regions[i, "start"]
  }

  return(mask_groups)
}

#' Compute Insertion CIGAR String
#'
#' @description Generate the CIGAR string representation for deletions.
#' @param position Position of deletion.
#' @param ref_position Position in the reference sequence.
#' @param ref_char Reference characters.
#' @param insertion_group Group identifier for deletion.
#' @return A string representing the CIGAR format for deletion.
compute_insertion_cigar <- function(position, ref_position, query_char, insertion_group) {
  # Determine the unique group value
  unique_group <- unique(insertion_group)[!is.na(unique(insertion_group))]

  # If there's no unique group value, return an NA result (this shouldn't happen, but just to be safe)
  if (length(unique_group) == 0) {
    return(NA_character_)
  }

  # Filter positions and characters based on the unique group value
  group_positions <- position[insertion_group == unique_group]
  group_ref_positions <- ref_position[insertion_group == unique_group]
  group_chars <- query_char[insertion_group == unique_group]

  # Get the starting position of the insertion relative to the reference
  start_position_ref <- ifelse(first(group_ref_positions) == 0, 1, first(group_ref_positions) + 1)

  # Return the pseudo CIGAR representation for the insertion
  return(paste0(start_position_ref, "I=", paste0(group_chars, collapse = "")))
}

#' Compute Deletion CIGAR String
#'
#' @description Generate the CIGAR string representation for deletions.
#' @param position Position of deletion.
#' @param ref_position Position in the reference sequence.
#' @param ref_char Reference characters.
#' @param deletion_group Group identifier for deletion.
#' @return A string representing the CIGAR format for deletion.
compute_deletion_cigar <- function(position, ref_position, ref_char, deletion_group) {
  # Determine the unique group value
  unique_group <- unique(deletion_group)[!is.na(unique(deletion_group))]

  # If there's no unique group value, return an NA result
  if (length(unique_group) == 0) {
    return(NA_character_)
  }

  # Filter positions and characters based on the unique group value
  group_positions <- position[deletion_group == unique_group]
  group_ref_positions <- ref_position[deletion_group == unique_group]
  group_chars <- ref_char[deletion_group == unique_group]

  # Get the starting position of the deletion relative to the reference
  start_position_ref <- first(group_ref_positions)

  # Return the pseudo CIGAR representation for the deletion with characters
  return(paste0(start_position_ref, "D=", paste0(group_chars, collapse = "")))
}

#' Compute Mask CIGAR String
#'
#' @description Generate the CIGAR string representation for mask regions.
#' @param position Position of the mask.
#' @param ref_position Position in the reference sequence.
#' @param ref_chars Characters from the reference sequence.
#' @param mask_group Group identifier for the mask region.
#' @return A string representing the CIGAR format for masks.
compute_mask_cigar <- function(position, ref_position, ref_chars, mask_group) {

  # Determine the unique group value
  unique_group <- unique(mask_group)[!is.na(unique(mask_group))]

  # If there's no unique group value, return an NA result
  if (length(unique_group) == 0) {
    return(NA_character_)
  }

  # Filter positions and ref_positions based on the unique group value
  group_positions <- position[mask_group == unique_group]
  group_ref_positions <- ref_position[mask_group == unique_group]
  group_ref_chars <- ref_chars[mask_group == unique_group]

  # Get the starting position of the mask relative to the reference
  start_position_ref <- first(group_ref_positions)

  # Calculate the length of the mask group and subtract the number of '-' in the group
  length_of_mask <- length(group_positions) - sum(group_ref_chars == "-")

  # Return the pseudo CIGAR representation for the mask
  return(paste0(start_position_ref, "+", length_of_mask, "N"))
}

#' Compute Substitution Group
#' 
#' @description Identify the group of substitution for given query and reference characters.
#' @param query_chars A vector of characters representing the query.
#' @param ref_chars A vector of characters representing the reference.
#' @param ref_position Position in the reference sequence.
#' @param mask_group A vector of integers representing mask groups.
#' @return A vector of integers representing substitution groups.
compute_substitution_group <- function(query_chars, ref_chars, ref_position, mask_group) {
  # ifelse(ref_chars != query_chars & ref_chars != "-" & query_chars != "-", ref_position, NA)
  if_else(!is.na(mask_group) | ref_chars == query_chars | ref_chars == "-" | query_chars == "-", NA, ref_position)
}


#' Build PseudoCIGAR String
#'
#' @description Create a pseudoCIGAR string based on the reference and query sequences.
#' @param reference Reference sequence.
#' @param query Query sequence.
#' @return A string representing the pseudoCIGAR format.
build_pseudoCIGAR_string <- function(reference, query) {

  # Check if lengths of reference and query are the same
  if (nchar(reference) != nchar(query)) {
    stop("The lengths of reference and query sequences must be the same.")
  }

  # Transform strings into vectors of characters
  ref_chars <- str_split(reference, pattern = "")[[1]]
  query_chars <- str_split(query, pattern = "")[[1]]

  # Check for positions where both sequences have a "-"
  dual_gaps <- ref_chars == "-" & query_chars == "-"
  if (any(dual_gaps)) {
    warning("Both reference and query sequences have '-' at the same position(s). This may indicate an issue with your data.")
  }

  # Prepare a dataframe with position level annotations of whether
  # there is an insertion, a deletion, a substitution, or whether the position is
  # in a masked region. Keep track of reference positions ('ref_position')
  # because all position reporting should be relative to the reference,
  # and we will need to keep track of this due to insertions in the ASV.
  df <- tibble(position = seq_along(ref_chars),
               ref_position = cumsum(ifelse(ref_chars != "-", 1, 0)),
               ref_char = ref_chars,
               query_char = query_chars,
               mask_group = compute_mask_group(ref_chars),
               insertion_group = compute_insertion_group(ref_chars, mask_group),
               deletion_group = compute_deletion_group(query_chars, mask_group),
               substitution_group = compute_substitution_group(query_chars, ref_chars, ref_position, mask_group)) %>%
    filter(!(is.na(mask_group) & is.na(insertion_group) & is.na(deletion_group) & is.na(substitution_group)))

  # Filter out the masked positions for insertions and deletions
  df <- df %>%
    mutate(
      insertion_group = ifelse(!is.na(mask_group), NA, insertion_group),
      deletion_group = ifelse(!is.na(mask_group), NA, deletion_group),
      substitution_group = ifelse(!is.na(mask_group), NA, substitution_group)
    )

  # CIGAR computation
  df_cigar <- df %>%
    group_by(mask_group, insertion_group, deletion_group, substitution_group) %>% # 
    reframe(
      start_position = first(position),
      result = case_when(
        !is.na(first(mask_group)) ~ compute_mask_cigar(position, ref_position, ref_char, mask_group),  # Modified this line to add ref_char
        !is.na(first(insertion_group)) & first(query_char) != "-" ~ compute_insertion_cigar(position, ref_position, query_char, insertion_group),
        !is.na(first(deletion_group)) & first(ref_char) != "-" ~ compute_deletion_cigar(position, ref_position, ref_char, deletion_group),
        !is.na(first(substitution_group)) ~ paste0(first(ref_position), first(query_char)),
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(result) & result != "") %>%
    arrange(start_position)

  # Combine the results into a single string
  cigar_str <- str_c(df_cigar$result, collapse = "")
  if (str_length(cigar_str) == 0) {
    cigar_str <- "." # if the cigar string is empty, replace with a dot
  }

  return(cigar_str)
}


# Main code to create the pseudocigar string
library(argparse)

parser <- ArgumentParser(description='Create the pseudoCIGAR string using masked or unmasked ASVs')
parser$add_argument('--alignments', type="character", required=TRUE, help = "File containing aligned ASVs")
parser$add_argument("--ncores", type="integer", default=1, help="Number of cores to use for parallel processing")

args <- parser$parse_args()
print(args)

# Load in alignment data. This can be masked (with 'N's)or unmasked data.
df.aln <- read.csv(args$alignments, sep="\t", header=TRUE)

# Setup parallel backend if needed
doMC::registerDoMC(args$ncores)

# Create the pseudoCIGAR string for each alignment entry
pseudo_cigar <- foreach(ii = 1:nrow(df.aln), .combine='bind_rows', .packages=c("stringr")) %dopar% {
  row <- df.aln[ii, ]
  tibble(
    sampleID = row$sampleID,
    refid = row$refid,
    asv = row$asv,
    pseudo_cigar=build_pseudoCIGAR_string(
      row$refseq,
      row$hapseq
    )
  )
}

# Write out the data frame
write.table(pseudo_cigar,file="alignments.pseudocigar.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
