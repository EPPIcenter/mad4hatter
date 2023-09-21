#!/usr/bin/env Rscript
################################################################################
# Script for Filtering Errors in Demultiplexed Amplicon Data
#
# This script is designed to process multiple samples that have been 
# demultiplexed using cutadapt. The script applies the DADA2 error filtering 
# process to each demultiplexed amplicon file.
#
# The primary goal is to improve the quality of the reads by removing those that 
# are likely erroneous. This is achieved by setting various thresholds for expected
# number of errors and length. 
#
# The final output is a "filtered" directory for each sample, containing the 
# filtered demultiplexed amplicons. There will also be a single metadata file 
# that provides a summary of the number of reads that were input and the number 
# that remain after filtering.
################################################################################

library(argparse)
library(dada2)

# Argument Parsing and Setup
parser <- ArgumentParser(description='Filter Errors')
parser$add_argument('--trimmed-path', nargs='+', type="character", required=TRUE, help="Path to trimmed fastq files")
parser$add_argument('--minLen', type="numeric", required=TRUE, help="Minimum length")
parser$add_argument('--maxN', type="numeric", required=TRUE, help="Maximum N")
parser$add_argument('--rm-phix', type="logical", required=TRUE, help="Remove PhiX")
parser$add_argument('--compress', type="logical", required=TRUE, help="Compress")
parser$add_argument('--ncores', type="numeric", required=TRUE, help="Number of threads to use")
parser$add_argument('--matchIDs', type="logical", required=TRUE, help="Match IDs")
parser$add_argument('--verbose', action="store_true", help="Verbose")

parser$add_argument('--trimRight_R1', type="numeric", required=TRUE, help="Trim Right")
parser$add_argument('--trimLeft_R1', type="numeric", required=TRUE, help="Trim Left")
parser$add_argument('--truncQ_R1', type="numeric", required=TRUE, help="Truncate quality")
parser$add_argument('--maxEE_R1', type="numeric", required=TRUE, help="Maximum expected error")

parser$add_argument('--trimRight_R2', type="numeric", required=TRUE, help="Trim Right")
parser$add_argument('--trimLeft_R2', type="numeric", required=TRUE, help="Trim Left")
parser$add_argument('--truncQ_R2', type="numeric", required=TRUE, help="Truncate quality")
parser$add_argument('--maxEE_R2', type="numeric", required=TRUE, help="Maximum expected error")

args <- parser$parse_args()
print(args)

# Extract directory names to create corresponding output directories
dir_names <- basename(args$trimmed_path)
output_paths <- sprintf("filtered_%s", dir_names)

# Create output directories
lapply(output_paths, function(x) dir.create(x, recursive=TRUE, showWarnings=FALSE))

# List Files and create filter paths in a vectorized manner
fnFs_list <- lapply(args$trimmed_path, function(tp) sort(list.files(path=tp, pattern="_R1.fastq.gz", recursive=TRUE, full.names = TRUE)))
fnRs_list <- lapply(args$trimmed_path, function(tp) sort(list.files(path=tp, pattern="_R2.fastq.gz", recursive=TRUE, full.names = TRUE)))

sample_names_list <- lapply(fnFs_list, function(fnFs) sapply(strsplit(basename(fnFs), "_R1"), `[`, 1))

filtFs_list <- mapply(function(op, sn) paste0(op, "/", sn, "_F_filt.fastq.gz"), output_paths, sample_names_list, SIMPLIFY = FALSE)
filtRs_list <- mapply(function(op, sn) paste0(op, "/", sn, "_R_filt.fastq.gz"), output_paths, sample_names_list, SIMPLIFY = FALSE)

# Perform Filtering and save filter metadata
filter_metadata_list <- mapply(function(fnFs, fnRs, filtFs, filtRs) {
  filterAndTrim(
        fnFs, filtFs, fnRs, filtRs,
        maxN=args$maxN, maxEE=c(args$maxEE_R1, args$maxEE_R2), truncQ=c(args$truncQ_R1, args$truncQ_R2),
        rm.phix=args$rm_phix, compress=args$compress, multithread=args$ncores, 
        trimRight=c(args$trimRight_R1, args$trimRight_R2), trimLeft=c(args$trimLeft_R1, args$trimLeft_R2), minLen=args$minLen, matchIDs=args$matchIDs
    )
}, fnFs_list, fnRs_list, filtFs_list, filtRs_list, SIMPLIFY = FALSE)

# Save filtered data and filter metadata
# Adjusted to save as "<sample_baseName>filter_metadata.RDS"
mapply(function(metadata, output) {
    output_filename <- file.path(output, sprintf("%s_filter_metadata.RDS", basename(output)))
    saveRDS(metadata, file=output_filename)
}, filter_metadata_list, output_paths)

