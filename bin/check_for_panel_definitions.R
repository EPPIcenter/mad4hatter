#!/usr/bin/env Rscript

library(BSgenome)
library(Biostrings)
library(argparse)
library(tidyverse)

parser <- ArgumentParser(prog = "Reference Validation Utility Script", description = "Utility to check if the reference contains definitions for the amplicons in the panel")
parser$add_argument("--reference-fasta", type = "character", help = "Artificial Reference FASTA", required = TRUE)
parser$add_argument("--ampliconFILE", type = "character", required = TRUE)

args <- parser$parse_args()

print(args)

# Function to compare the declared panel and the definitions of those amplicons
# in the reference sequence
check_definitions_for_amplicon_panel <- function(reference_filepath, panel_filepath) {
  panel_df <- read_table(panel_filepath, show_col_types = FALSE)
  ref_sequences <- Biostrings::readDNAStringSet(reference_filepath)

  check_for_invalid_definitions <- function() {
    amplicon_in_reference <- names(ref_sequences)
    amplicon_width_in_reference <- width(ref_sequences)
    panel_amplicon_declarations <- panel_df$amplicon

    # Make sure the amplicon is defined in the reference - definitions _cannot_ be empty.
    amplicon_in_reference[!amplicon_in_reference %in% panel_amplicon_declarations | amplicon_width_in_reference == 0]
  }

  results <- check_for_invalid_definitions()
  
  return (results)
}

missing_amplicons <- check_definitions_for_amplicon_panel(args$reference_fasta, args$ampliconFILE)

any_missing_amplicons <- function(results) {
  length(missing_amplicons) > 0
}

# Declare error codes to return at end of script
ERROR_SUCCESS <- 0
ERROR_FAILURE <- 1

if (any_missing_amplicons(missing_amplicons)) {
  errmsg <- paste0("Amplicons with missing definitions detected: [", paste(missing_amplicons, collapse = ", "), "].")

  # Output to stderr will cause Nextflow to shutdown.
  write(errmsg, stderr())

  # Additionally, return an error code.
  errcode <- ERROR_FAILURE

} else {
  print("No missing amplicons detected in reference!")
  errcode <- ERROR_SUCCESS
}

# Return error code
quit(status = errcode)
