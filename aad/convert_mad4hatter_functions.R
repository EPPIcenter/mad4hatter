library(tidyverse)

# -----------------------------------------------------------------------------
# release_version DETECTION
# -----------------------------------------------------------------------------

detect_release_version <- function(input_dir) {
  
  files <- list.files(input_dir, recursive = TRUE)
  
  # Read allele_data.txt headers to fingerprint release_version
  allele_file <- file.path(input_dir, "allele_data.txt")
  if (!file.exists(allele_file)) {
    stop("Cannot find allele_data.txt in input_dir. Is this a mad4hatter results folder?")
  }
  
  allele_cols <- names(read_tsv(allele_file, n_max = 0, show_col_types = FALSE))
  
  # v0.1.8: lowercase sampleID, lowercase locus, lowercase asv, has 'allele' column
  if ("sampleID" %in% allele_cols && "locus" %in% allele_cols && "allele" %in% allele_cols) {
    release_version <- "0.1.8"
    
    # v0.2.2: uppercase SampleID, uppercase Locus, has PseudoCIGAR (no 'allele')  
  } else if ("SampleID" %in% allele_cols && "Locus" %in% allele_cols && "PseudoCIGAR" %in% allele_cols) {
    release_version <- "0.2.2"
    
    # v1.0.0: sample_name, target_name, pseudocigar_unmasked
  } else if ("sample_name" %in% allele_cols && "target_name" %in% allele_cols) {
    release_version <- "1.0.0"
    
  } else {
    stop(paste(
      "Could not detect release_version from allele_data.txt column names.",
      "Found columns:", paste(allele_cols, collapse = ", ")
    ))
  }
  
  message(sprintf("Detected release_version: %s", release_version))
  return(release_version)
}

# -----------------------------------------------------------------------------
# LOCUS LOOKUP HELPER
# -----------------------------------------------------------------------------

apply_locus_lookup <- function(df, locus_col, lookup) {
  # dplyr::rename the locus column to old_name for joining
  df <- df %>% dplyr::rename(old_name = !!sym(locus_col))
  
  # Check for unmatched loci
  unmatched <- setdiff(unique(df$old_name), lookup$old_name)
  if (length(unmatched) > 0) {
    warning(sprintf(
      "The following loci were NOT found in the lookup file and will be kept as-is with pool = NA:\n  %s",
      paste(unmatched, collapse = "\n  ")
    ))
  }
  
  df <- df %>%
    left_join(lookup, by = "old_name") %>%
    mutate(
      target_name = if_else(is.na(target_name), old_name, target_name),
      pool        = if_else(old_name %in% unmatched, NA_character_, pool)
    ) %>%
    select(-old_name)
  
  return(df)
}

# -----------------------------------------------------------------------------
# FILE CONVERTERS
# -----------------------------------------------------------------------------

convert_allele_data <- function(input_dir, output_dir, release_version, lookup) {
  f <- file.path(input_dir, "allele_data.txt")
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (release_version == "0.1.8") {
    df <- df %>%
      dplyr::rename(
        sample_name      = sampleID,
        asv              = asv,
        reads            = reads,
        pseudocigar_masked = pseudo_cigar
      ) %>%
      select(-allele) %>%  # drop allele column
      mutate(asv_masked           = NA_character_,
             pseudocigar_unmasked = NA_character_)
    locus_col <- "locus"
    
  } else if (release_version == "0.2.2") {
    df <- df %>%
      dplyr::rename(
        sample_name        = SampleID,
        asv                = ASV,
        reads              = Reads,
        pseudocigar_masked = PseudoCIGAR
      ) %>%
      mutate(asv_masked           = NA_character_,
             pseudocigar_unmasked = NA_character_)
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
  df <- allele_data %>%
    group_by(sample_name, target_name, asv_masked, pseudocigar_masked, pool) %>%
    summarise(reads = sum(reads, na.rm = TRUE), .groups = "drop")
  
  # Final column order matching v1.0.0
  df <- df %>%
    select(sample_name, target_name,
           asv_masked, pseudocigar_masked, reads, pool)
  
  write_tsv(df, file.path(output_dir, "allele_data_collapsed.txt"))
  message("  [OK] allele_data_collapsed.txt")
}

convert_sample_coverage <- function(input_dir, output_dir, release_version) {
  f <- file.path(input_dir, "sample_coverage.txt")
  df <- read_tsv(f, show_col_types = FALSE)
  
  # Both v0.1.8 and v0.2.2 have same columns, just need sample_name dplyr::rename
  df <- df %>%
    dplyr::rename(sample_name = SampleID,
                  stage       = Stage,
                  reads       = Reads)
  
  write_tsv(df, file.path(output_dir, "sample_coverage.txt"))
  message("  [OK] sample_coverage.txt")
}

convert_amplicon_coverage <- function(input_dir, output_dir, release_version, lookup) {
  f <- file.path(input_dir, "amplicon_coverage.txt")
  df <- read_tsv(f, show_col_types = FALSE)
  
  # Both old release_versions have same column names
  df <- df %>%
    dplyr::rename(sample_name = SampleID,
                  reads       = Reads)
  
  # Apply locus lookup
  df <- apply_locus_lookup(df, "Locus", lookup) %>%
    select(sample_name, target_name, reads, OutputDada2, OutputPostprocessing)
  
  write_tsv(df, file.path(output_dir, "amplicon_coverage.txt"))
  message("  [OK] amplicon_coverage.txt")
}

convert_resmarker_table <- function(input_dir, output_dir, release_version, lookup) {
  f <- file.path(input_dir, "resistance_marker_module", "resmarker_table.txt")
  if (!file.exists(f)) return(invisible(NULL))
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (release_version == "0.1.8") {
    df <- df %>%
      dplyr::rename(
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
      group_by(sample_name, gene_id, gene, aa_position, ref_codon, codon, codon_ref_alt, ref_aa, aa, aa_ref_alt) %>%  # CodonoStart dropped in later release_versions
      summarize(reads = sum(reads,na.rm = F),
                .groups = "drop") %>% 
      mutate(follows_indel  = NA,
             codon_masked   = NA,
             multiple_loci  = NA)
    
  } else if (release_version == "0.2.2") {
    df <- df %>%
      dplyr::rename(
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

convert_resmarker_table_by_locus <- function(input_dir, output_dir, release_version, lookup) {
  f <- file.path(input_dir, "resistance_marker_module", "resmarker_table_by_locus.txt")
  if (!file.exists(f)) {
    message(
      "  [SKIP] resistance_marker_module/resmarker_table_by_locus.txt does not exist - no output will be created"
    )
    return(invisible(NULL))
  }
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (release_version == "0.1.8") {
    # This file doesn't exist in v0.1.8, skip
    return(invisible(NULL))
  } else if (release_version == "0.2.2") {
    df <- df %>%
      dplyr::rename(
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

convert_resmarker_microhaplotype <- function(input_dir, output_dir, release_version, lookup) {
  
  # File name differs by release_version
  if (release_version == "0.1.8") {
    f <- file.path(input_dir, "resistance_marker_module", "resmarker_microhap_table.txt")
  } else {
    f <- file.path(input_dir, "resistance_marker_module", "resmarker_microhaplotype_table.txt")
  }
  if (!file.exists(f)) return(invisible(NULL))
  df <- read_tsv(f, show_col_types = FALSE)
  
  if (release_version == "0.1.8") {
    # MicrohapIndex -> mhap_aa_positions
    # RefMicrohap   -> ref_mhap
    # Microhaplotype -> mhap  (note: MicrohapRefAlt and Reads are separate columns)
    df <- df %>%
      dplyr::rename(
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
    
  } else if (release_version == "0.2.2") {
    df <- df %>%
      dplyr::rename(
        sample_name       = SampleID,
        gene_id           = GeneID,
        gene              = Gene,
        mhap_aa_positions = MicrohaplotypeCodonIDs,
        ref_mhap          = RefMicrohap,
        mhap              = Microhaplotype,
        mhap_ref_alt      = MicrohapRefAlt,
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

convert_all_mutations <- function(input_dir, output_dir, release_version, lookup) {
  # Per spec: resmarker_new_mutations.txt (v0.1.8) and all_mutations_table.txt (v0.2.2)
  # should both be DELETED (not converted) in the output
  
  old_files <- c(
    file.path(input_dir, "resistance_marker_module", "resmarker_new_mutations.txt"),
    file.path(input_dir, "resistance_marker_module", "all_mutations_table.txt")
  )
  for (f in old_files) {
    if (file.exists(f)) {
      message(sprintf("  [SKIP] %s — this file is intentionally excluded as versions prior to v1.0.0 were not correct", basename(f)))
    }
  }
}

copy_passthrough_files <- function(input_dir, output_dir, release_version) {
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

create_panel_information <- function(output_dir, amplicon_info, references, resmarker_info) {
  
  panel_dir <- file.path(output_dir, "panel_information")
  dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get target_names present in this run from the already-converted amplicon_coverage.txt
  targets_in_data <- read_tsv(file.path(output_dir, "amplicon_coverage.txt"),
                              show_col_types = FALSE) %>%
    pull(target_name) %>%
    unique()
  
  # --- amplicon_info.tsv -----------------------------------------------------
  if (!file.exists(amplicon_info)) {
    warning("amplicon_info file not found, skipping: ", amplicon_info)
  } else {
    df_amplicon <- read_tsv(amplicon_info, show_col_types = FALSE)
    
    missing_targets <- setdiff(targets_in_data, df_amplicon$target_name)
    if (length(missing_targets) > 0) {
      warning(sprintf(
        "The following target_names from the data were NOT found in amplicon_info:\n  %s",
        paste(missing_targets, collapse = "\n  ")
      ))
    }
    
    df_amplicon %>%
      filter(target_name %in% targets_in_data) %>%
      write_tsv(file.path(panel_dir, "amplicon_info.tsv"))
    message("  [OK] panel_information/amplicon_info.tsv")
  }
  
  # --- reference.fasta -------------------------------------------------------
  if (!file.exists(references)) {
    warning("references fasta file not found, skipping: ", references)
  } else {
    fasta_lines <- readLines(references)
    
    # Find all header lines and extract sequence names (everything after > up to first space)
    header_idx <- which(startsWith(fasta_lines, ">"))
    seq_names  <- sub("^>([^ ]+).*", "\\1", fasta_lines[header_idx])
    
    # Build per-entry line ranges
    entry_start <- header_idx
    entry_end   <- c(header_idx[-1] - 1, length(fasta_lines))
    
    missing_refs <- setdiff(targets_in_data, seq_names)
    if (length(missing_refs) > 0) {
      warning(sprintf(
        "The following target_names from the data were NOT found in the reference fasta:\n  %s",
        paste(missing_refs, collapse = "\n  ")
      ))
    }
    
    # Keep only entries whose name matches a target in the data
    keep      <- seq_names %in% targets_in_data
    kept_lines <- unlist(mapply(
      function(s, e) fasta_lines[s:e],
      entry_start[keep], entry_end[keep],
      SIMPLIFY = FALSE
    ))
    
    writeLines(kept_lines, file.path(panel_dir, "reference.fasta"))
    message("  [OK] panel_information/reference.fasta")
  }
  
  # --- resmarker_info.tsv ----------------------------------------------------
  if (!file.exists(resmarker_info)) {
    warning("resmarker_info file not found, skipping: ", resmarker_info)
  } else {
    resmarker_table_path <- file.path(output_dir, "resistance_marker_module", "resmarker_table.txt")
    
    if (!file.exists(resmarker_table_path)) {
      message("  [SKIP] panel_information/resmarker_info.tsv — no resmarker_table.txt in output")
    } else {
      # Keys present in this run's resmarker output
      resmarker_keys <- read_tsv(resmarker_table_path, show_col_types = FALSE) %>%
        distinct(gene_id, gene, aa_position)
      
      df_resmarker <- read_tsv(resmarker_info, show_col_types = FALSE)
      
      # Warn about any keys in the data that are absent from the reference file
      missing_keys <- anti_join(resmarker_keys, df_resmarker,
                                by = c("gene_id", "gene", "aa_position"))
      if (nrow(missing_keys) > 0) {
        warning(sprintf(
          "The following gene_id/gene/aa_position combinations from resmarker_table.txt were NOT found in resmarker_info:\n  %s",
          paste(apply(missing_keys, 1, paste, collapse = " / "), collapse = "\n  ")
        ))
      }
      
      df_resmarker %>%
        semi_join(resmarker_keys, by = c("gene_id", "gene", "aa_position")) %>%
        mutate(target_name = NA_character_,
               codon_start_in_target = NA_character_) %>% 
        write_tsv(file.path(panel_dir, "resmarker_info.tsv"))
      message("  [OK] panel_information/resmarker_info.tsv")
    }
  }
}



# -----------------------------------------------------------------------------
# MAIN FUNCTION
# -----------------------------------------------------------------------------
convert_mad4hatter <- function(input_dir, output_dir, locus_lookup,
                               amplicon_info, references, resmarker_info,
                               release_version = NULL) {
  # --- Validate inputs -------------------------------------------------------
  if (!dir.exists(input_dir))  stop("input_dir does not exist: ", input_dir)
  if (!file.exists(locus_lookup)) stop("locus_lookup file not found: ", locus_lookup)
  
  # --- Detect release_version --------------------------------------------------------
  if (is.null(release_version)) {
    release_version <- detect_release_version(input_dir)
  } else {
    message(sprintf("Using user-supplied release_version: %s", release_version))
  }
  
  if (release_version == "1.0.0") {
    message("Input already appears to be v1.0.0. No conrelease_version needed.")
    return(invisible(NULL))
  }
  
  if (!release_version %in% c("0.1.8", "0.2.2")) {
    stop(sprintf("Unsupported release_version '%s'. Supported: 0.1.8, 0.2.2", release_version))
  }
  
  # --- Load lookup -----------------------------------------------------------
  lookup <- read_tsv(locus_lookup, show_col_types = FALSE)
  required_cols <- c("old_name", "new_name", "pool")
  missing_cols  <- setdiff(required_cols, names(lookup))
  if (length(missing_cols) > 0) {
    stop("locus_lookup is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  lookup <- lookup %>% 
    dplyr::rename(target_name = new_name) 
  
  # --- Create output directories ---------------------------------------------
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "resistance_marker_module"), showWarnings = FALSE)
  
  message(sprintf("\nConverting from v%s → v1.0.0", release_version))
  message(sprintf("Input:  %s", input_dir))
  message(sprintf("Output: %s\n", output_dir))
  
  # --- Convert files ---------------------------------------------------------
  allele_data <- convert_allele_data(input_dir, output_dir, release_version, lookup)
  convert_allele_data_collapsed(allele_data, output_dir)
  convert_sample_coverage(input_dir, output_dir, release_version)
  convert_amplicon_coverage(input_dir, output_dir, release_version, lookup)
  convert_resmarker_table(input_dir, output_dir, release_version, lookup)
  convert_resmarker_table_by_locus(input_dir, output_dir, release_version, lookup)
  convert_resmarker_microhaplotype(input_dir, output_dir, release_version, lookup)
  convert_all_mutations(input_dir, output_dir, release_version, lookup)
  copy_passthrough_files(input_dir, output_dir, release_version)
  
  # --- Create panel_information ----------------------------------------------
  create_panel_information(output_dir, amplicon_info, references, resmarker_info)
  
  
  message(sprintf("\nDone. Output written to: %s", output_dir))
  invisible(output_dir)
}
