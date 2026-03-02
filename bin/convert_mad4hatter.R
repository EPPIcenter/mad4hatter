# =============================================================================
# convert_mad4hatter.R
# Converts mad4hatter results from v0.1.8 or v0.2.2 to v1.0.0 format
#
# Usage:
#   source("convert_mad4hatter.R")
#   convert_mad4hatter(
#     input_dir   = "path/to/old/results",
#     output_dir  = "path/to/new/results",
#     locus_lookup = "path/to/locus_lookup.csv"
#   )
#
# locus_lookup.csv must have columns: old_target_name, target_name, pool
# =============================================================================

library(dplyr)
library(readr)
library(stringr)

# -----------------------------------------------------------------------------
# VERSION DETECTION
# -----------------------------------------------------------------------------

detect_version <- function(input_dir) {
  
  files <- list.files(input_dir, recursive = TRUE)
  
  # Read allele_data.txt headers to fingerprint version
  allele_file <- file.path(input_dir, "allele_data.txt")
  if (!file.exists(allele_file)) {
    stop("Cannot find allele_data.txt in input_dir. Is this a mad4hatter results folder?")
  }
  
  allele_cols <- names(read_tsv(allele_file, n_max = 0, show_col_types = FALSE))
  
  # v0.1.8: lowercase sampleID, lowercase locus, lowercase asv, has 'allele' column
  if ("sampleID" %in% allele_cols && "locus" %in% allele_cols && "allele" %in% allele_cols) {
    version <- "0.1.8"
    
  # v0.2.2: uppercase SampleID, uppercase Locus, has PseudoCIGAR (no 'allele')  
  } else if ("SampleID" %in% allele_cols && "Locus" %in% allele_cols && "PseudoCIGAR" %in% allele_cols) {
    version <- "0.2.2"
    
  # v1.0.0: sample_name, target_name, pseudocigar_unmasked
  } else if ("sample_name" %in% allele_cols && "target_name" %in% allele_cols) {
    version <- "1.0.0"
    
  } else {
    stop(paste(
      "Could not detect version from allele_data.txt column names.",
      "Found columns:", paste(allele_cols, collapse = ", ")
    ))
  }
  
  message(sprintf("Detected version: %s", version))
  return(version)
}

# -----------------------------------------------------------------------------
# LOCUS LOOKUP HELPER
# -----------------------------------------------------------------------------

apply_locus_lookup <- function(df, locus_col, lookup) {
  # Rename the locus column to old_target_name for joining
  df <- df %>% rename(old_target_name = !!sym(locus_col))
  
  # Check for unmatched loci
  unmatched <- setdiff(unique(df$old_target_name), lookup$old_target_name)
  if (length(unmatched) > 0) {
    warning(sprintf(
      "The following loci were NOT found in the lookup file and will be kept as-is with pool = NA:\n  %s",
      paste(unmatched, collapse = "\n  ")
    ))
  }
  
  df <- df %>%
    left_join(lookup, by = "old_target_name") %>%
    mutate(
      target_name = if_else(is.na(target_name), old_target_name, target_name),
      pool        = if_else(old_target_name %in% unmatched, NA_character_, pool)
    ) %>%
    select(-old_target_name)
  
  return(df)
}

# -----------------------------------------------------------------------------
# FILE CONVERTERS
# -----------------------------------------------------------------------------

convert_allele_data <- function(input_dir, output_dir, version, lookup) {
  f <- file.path(input_dir, "allele_data.txt")
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (version == "0.1.8") {
    df <- df %>%
      rename(
        sample_name      = sampleID,
        asv              = asv,
        reads            = reads,
        pseudocigar_masked = pseudo_cigar
      ) %>%
      select(-allele) %>%  # drop allele column
      mutate(asv_masked           = NA_character_,
             pseudocigar_unmasked = "FIX_THIS")
    locus_col <- "locus"
    
  } else if (version == "0.2.2") {
    df <- df %>%
      rename(
        sample_name        = SampleID,
        asv                = ASV,
        reads              = Reads,
        pseudocigar_masked = PseudoCIGAR
      ) %>%
      mutate(asv_masked           = NA_character_,
             pseudocigar_unmasked = "FIX_THIS")
    locus_col <- "Locus"
  }
  
  # Apply locus lookup (adds target_name and pool, removes locus column)
  df <- apply_locus_lookup(df, locus_col, lookup)
  
  # Final column order matching v1.0.0
  df <- df %>%
    select(sample_name, target_name, asv, pseudocigar_unmasked,
           asv_masked, pseudocigar_masked, reads, pool)
  
  write_tsv(df, file.path(output_dir, "allele_data.txt"))
  message("  [OK] allele_data.txt")
  return(df)
}

convert_allele_data_collapsed <- function(allele_data, output_dir) {
  # Derived from allele_data: group by sample/target/asv_masked/pseudocigar_masked/pool
  # Remove rows where asv_masked is NA (can't collapse without it)
  df <- allele_data %>%
    filter(!is.na(asv_masked), !is.na(pseudocigar_masked)) %>%
    group_by(sample_name, target_name, asv_masked, pseudocigar_masked, pool) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop")
  
  write_tsv(df, file.path(output_dir, "allele_data_collapsed.txt"))
  message("  [OK] allele_data_collapsed.txt (note: rows with NA asv_masked were excluded)")
}

convert_sample_coverage <- function(input_dir, output_dir, version) {
  f <- file.path(input_dir, "sample_coverage.txt")
  df <- read_tsv(f, show_col_types = FALSE)
  
  # Both v0.1.8 and v0.2.2 have same columns, just need sample_name rename
  df <- df %>%
    rename(sample_name = SampleID,
           stage       = Stage,
           reads       = Reads)
  
  write_tsv(df, file.path(output_dir, "sample_coverage.txt"))
  message("  [OK] sample_coverage.txt")
}

convert_amplicon_coverage <- function(input_dir, output_dir, version, lookup) {
  f <- file.path(input_dir, "amplicon_coverage.txt")
  df <- read_tsv(f, show_col_types = FALSE)
  
  # Both old versions have same column names
  df <- df %>%
    rename(sample_name = SampleID,
           reads       = Reads)
  
  # Apply locus lookup
  df <- apply_locus_lookup(df, "Locus", lookup) %>%
    select(sample_name, target_name, reads, OutputDada2, OutputPostprocessing)
  
  write_tsv(df, file.path(output_dir, "amplicon_coverage.txt"))
  message("  [OK] amplicon_coverage.txt")
}

convert_resmarker_table <- function(input_dir, output_dir, version, lookup) {
  f <- file.path(input_dir, "resistance_marker_module", "resmarker_table.txt")
  if (!file.exists(f)) return(invisible(NULL))
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (version == "0.1.8") {
    df <- df %>%
      rename(
        sample_name  = SampleID,
        gene_id      = GeneID,
        gene         = Gene,
        aa_position  = CodonID,
        ref_codon    = RefCodon,
        codon        = Codon,
        codon_ref_alt = CodonRefAlt,
        ref_aa       = RefAA,
        aa           = AA,
        aa_ref_alt   = AARefAlt,
        reads        = Reads
      ) %>%
      select(-CodonStart) %>%  # dropped in later versions
      mutate(follows_indel  = NA,
             codon_masked   = NA,
             multiple_loci  = NA)
    
  } else if (version == "0.2.2") {
    df <- df %>%
      rename(
        sample_name   = SampleID,
        gene_id       = GeneID,
        gene          = Gene,
        aa_position   = CodonID,
        ref_codon     = RefCodon,
        codon         = Codon,
        codon_ref_alt = CodonRefAlt,
        ref_aa        = RefAA,
        aa            = AA,
        aa_ref_alt    = AARefAlt,
        follows_indel = FollowsIndel,
        codon_masked  = CodonMasked,
        multiple_loci = MultipleLoci,
        reads         = Reads
      )
  }
  
  df <- df %>%
    select(sample_name, gene_id, gene, aa_position, ref_codon, codon,
           codon_ref_alt, ref_aa, aa, aa_ref_alt, follows_indel,
           codon_masked, multiple_loci, reads)
  
  write_tsv(df, file.path(output_dir, "resistance_marker_module", "resmarker_table.txt"))
  message("  [OK] resistance_marker_module/resmarker_table.txt")
}

convert_resmarker_table_by_locus <- function(input_dir, output_dir, version, lookup) {
  f <- file.path(input_dir, "resistance_marker_module", "resmarker_table_by_locus.txt")
  if (!file.exists(f)) return(invisible(NULL))
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (version == "0.1.8") {
    # This file doesn't exist in v0.1.8, skip
    return(invisible(NULL))
  } else if (version == "0.2.2") {
    df <- df %>%
      rename(
        sample_name   = SampleID,
        gene_id       = GeneID,
        gene          = Gene,
        aa_position   = CodonID,
        ref_codon     = RefCodon,
        codon         = Codon,
        codon_ref_alt = CodonRefAlt,
        ref_aa        = RefAAAA,   # note: typo in v0.2.2 header
        aa            = AA,
        aa_ref_alt    = AARefAlt,
        follows_indel = FollowsIndel,
        codon_masked  = CodonMasked,
        reads         = Reads
      )
    df <- apply_locus_lookup(df, "Locus", lookup)
  }
  
  df <- df %>%
    select(sample_name, gene_id, gene, target_name, aa_position, ref_codon,
           codon, codon_ref_alt, ref_aa, aa, aa_ref_alt, follows_indel,
           codon_masked, reads)
  
  write_tsv(df, file.path(output_dir, "resistance_marker_module", "resmarker_table_by_locus.txt"))
  message("  [OK] resistance_marker_module/resmarker_table_by_locus.txt")
}

convert_resmarker_microhaplotype <- function(input_dir, output_dir, version, lookup) {
  
  # File name differs by version
  if (version == "0.1.8") {
    f <- file.path(input_dir, "resistance_marker_module", "resmarker_microhap_table.txt")
  } else {
    f <- file.path(input_dir, "resistance_marker_module", "resmarker_microhaplotype_table.txt")
  }
  if (!file.exists(f)) return(invisible(NULL))
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (version == "0.1.8") {
    # MicrohapIndex -> mhap_aa_positions
    # RefMicrohap   -> ref_mhap
    # Microhaplotype -> mhap  (note: MicrohapRefAlt and Reads are separate columns)
    df <- df %>%
      rename(
        sample_name      = SampleID,
        gene_id          = GeneID,
        gene             = Gene,
        mhap_aa_positions = MicrohapIndex,
        ref_mhap         = RefMicrohap,
        mhap             = Microhaplotype,
        mhap_ref_alt     = MicrohapRefAlt,
        reads            = Reads
      )
    # No Locus column in v0.1.8 microhap table, so target_name will be NA
    df <- df %>% mutate(target_name = NA_character_)
    
  } else if (version == "0.2.2") {
    df <- df %>%
      rename(
        sample_name       = SampleID,
        gene_id           = GeneID,
        gene              = Gene,
        mhap_aa_positions = MicrohaplotypeCodonIDs,
        ref_mhap          = RefMicrohap,
        mhap              = Microhaplotype,
        mhap_ref_alt      = MicrohaplotypeMicrohapRefAlt,
        reads             = Reads
      )
    df <- apply_locus_lookup(df, "Locus", lookup)
  }
  
  df <- df %>%
    select(sample_name, gene_id, gene, target_name, mhap_aa_positions,
           ref_mhap, mhap, mhap_ref_alt, reads)
  
  write_tsv(df, file.path(output_dir, "resistance_marker_module", "resmarker_microhaplotype_table.txt"))
  message("  [OK] resistance_marker_module/resmarker_microhaplotype_table.txt")
}

convert_all_mutations <- function(input_dir, output_dir, version, lookup) {
  # Per spec: resmarker_new_mutations.txt (v0.1.8) and all_mutations_table.txt (v0.2.2)
  # should both be DELETED (not converted) in the output
  
  old_files <- c(
    file.path(input_dir, "resistance_marker_module", "resmarker_new_mutations.txt"),
    file.path(input_dir, "resistance_marker_module", "all_mutations_table.txt")
  )
  for (f in old_files) {
    if (file.exists(f)) {
      message(sprintf("  [SKIP] %s — this file is intentionally excluded from v1.0.0 output", basename(f)))
    }
  }
}

copy_passthrough_files <- function(input_dir, output_dir, version) {
  # Files to copy as-is (quality_report, run, raw_dada2_output if present)
  passthrough_dirs <- c("quality_report", "run", "raw_dada2_output")
  
  for (d in passthrough_dirs) {
    src <- file.path(input_dir, d)
    dst <- file.path(output_dir, d)
    if (dir.exists(src)) {
      dir.create(dst, showWarnings = FALSE, recursive = TRUE)
      files <- list.files(src, full.names = TRUE)
      file.copy(files, dst)
      message(sprintf("  [COPY] %s/ (passed through unchanged)", d))
    }
  }
}

# -----------------------------------------------------------------------------
# MAIN FUNCTION
# -----------------------------------------------------------------------------

#' Convert mad4hatter results to v1.0.0 format
#'
#' @param input_dir   Path to the old results folder
#' @param output_dir  Path where v1.0.0-formatted results will be written
#' @param locus_lookup Path to CSV with columns: old_target_name, target_name, pool
#' @param version     Optional. If NULL (default), version is auto-detected.
#'                    Override with "0.1.8" or "0.2.2" if needed.
#'
#' @examples
#' convert_mad4hatter(
#'   input_dir    = "~/results/MULE_run1",
#'   output_dir   = "~/results/MULE_run1_v1",
#'   locus_lookup = "~/locus_lookup.csv"
#' )
convert_mad4hatter <- function(input_dir, output_dir, locus_lookup, version = NULL) {
  
  # --- Validate inputs -------------------------------------------------------
  if (!dir.exists(input_dir))  stop("input_dir does not exist: ", input_dir)
  if (!file.exists(locus_lookup)) stop("locus_lookup file not found: ", locus_lookup)
  
  # --- Detect version --------------------------------------------------------
  if (is.null(version)) {
    version <- detect_version(input_dir)
  } else {
    message(sprintf("Using user-supplied version: %s", version))
  }
  
  if (version == "1.0.0") {
    message("Input already appears to be v1.0.0. No conversion needed.")
    return(invisible(NULL))
  }
  
  if (!version %in% c("0.1.8", "0.2.2")) {
    stop(sprintf("Unsupported version '%s'. Supported: 0.1.8, 0.2.2", version))
  }
  
  # --- Load lookup -----------------------------------------------------------
  lookup <- read_csv(locus_lookup, show_col_types = FALSE)
  required_cols <- c("old_target_name", "target_name", "pool")
  missing_cols  <- setdiff(required_cols, names(lookup))
  if (length(missing_cols) > 0) {
    stop("locus_lookup is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # --- Create output directories ---------------------------------------------
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "resistance_marker_module"), showWarnings = FALSE)
  
  message(sprintf("\nConverting from v%s → v1.0.0", version))
  message(sprintf("  Input:  %s", input_dir))
  message(sprintf("  Output: %s\n", output_dir))
  
  # --- Convert files ---------------------------------------------------------
  allele_data <- convert_allele_data(input_dir, output_dir, version, lookup)
  convert_allele_data_collapsed(allele_data, output_dir)
  convert_sample_coverage(input_dir, output_dir, version)
  convert_amplicon_coverage(input_dir, output_dir, version, lookup)
  convert_resmarker_table(input_dir, output_dir, version, lookup)
  convert_resmarker_table_by_locus(input_dir, output_dir, version, lookup)
  convert_resmarker_microhaplotype(input_dir, output_dir, version, lookup)
  convert_all_mutations(input_dir, output_dir, version, lookup)
  copy_passthrough_files(input_dir, output_dir, version)
  
  message(sprintf("\nDone. Output written to: %s", output_dir))
  message("NOTE: 'pseudocigar_unmasked' columns are filled with 'FIX_THIS' and require manual computation.")
  
  invisible(output_dir)
}
