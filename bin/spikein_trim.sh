#!/usr/bin/env bash

# Function to display usage help
usage() {
    echo "Usage: $0 -f FASTQ_DIR -p PRIMER_CSV -1 FWD_PRIMER -2 REV_PRIMER -t THREADS (default: 4)"
    exit 1
}

# FASTQ_DIR should be the directory containing subdirectories of fastq.gz files. This
# allows for the creation of bam files from different trimming stages if desired.

# Function to handle SIGINT (Ctrl+C)
cleanup_and_exit() {
    echo "Script interrupted. Exiting..."
    exit 1
}

# Trap SIGINT and call cleanup_and_exit function
trap cleanup_and_exit SIGINT

# Default threads to use with cutadapt
THREADS=4

# Parse command-line options
while getopts "f:1:2:r:t:" opt; do
    case $opt in
        f) FASTQ_DIR=$OPTARG ;;
        1) FWD_PRIMER=$OPTARG ;;
        2) REV_PRIMER=$OPTARG ;;
        r) RESULTS_DIR=$OPTARG ;;
        t) THREADS=$OPTARG ;;
        *) usage ;;
    esac
done

# Check if required variables are set
if [ -z "$FASTQ_DIR" ] || [ -z "$RESULTS_DIR" ]; then
    echo "FASTQ_DIR and RESULTS_DIR must be set."
    usage
fi

# Function to trim adapters and demultiplex using cutadapt
trim_and_demultiplex_paired() {
    local fwd_fastq_file=$1
    local rev_fastq_file=$2
    local base_name=$(basename "$fwd_fastq_file" _unknown_R1.fastq.gz)
    local unknown_dir="$RESULTS_DIR/unknown/${base_name}"
    local demultiplexed_dir="$RESULTS_DIR/demultiplexed/${base_name}"
    local primer_json_report="$RESULTS_DIR/cutadapt_json/${base_name}_primers.json"
    local allowed_errors=0

    # Create output directories
    mkdir -p "$unknown_dir"
    mkdir -p "$demultiplexed_dir"
    mkdir -p "$RESULTS_DIR/cutadapt_json"

    # Check if FASTA files exist
    if [ ! -f "$FWD_PRIMER" ] || [ ! -f "$REV_PRIMER" ]; then
        echo "FASTA files not found. Exiting."
        exit 1
    fi

    # Run cutadapt for primer trimming, demultiplexing, and JSON reporting
    cutadapt -g file:"$FWD_PRIMER" \
             -G file:"$REV_PRIMER" \
             --trim-n -q 10 \
             --minimum-length 100 \
             -o "${demultiplexed_dir}/{name}_${base_name}_R1.fastq.gz" \
             -p "${demultiplexed_dir}/{name}_${base_name}_R2.fastq.gz" \
             -e ${allowed_errors} \
             --no-indels \
             --untrimmed-output ${unknown_dir}/${base_name}_unknown_R1.fastq.gz \
             --untrimmed-paired-output ${unknown_dir}/${base_name}_unknown_R2.fastq.gz \
             -j $THREADS \
             --json "$primer_json_report" \
             --quiet \
             "$fwd_fastq_file" \
             "$rev_fastq_file" > /dev/null
}

# Check if the FASTQ_DIR exists
if [ ! -d "$FASTQ_DIR" ]; then
    echo "FASTQ directory does not exist: $FASTQ_DIR"
    exit 1
fi

# Modify the loop for processing FASTQ files
for fwd_fastq_file in "$FASTQ_DIR"/*unknown_R1.fastq.gz; do
    # Construct the path to the corresponding reverse read file
    rev_fastq_file="${fwd_fastq_file/_R1/_R2}"

    # Check if the reverse read file exists
    if [ ! -f "$rev_fastq_file" ]; then
        echo "Reverse read file not found for $fwd_fastq_file"
        continue  # Skip to the next file
    fi

    trim_and_demultiplex_paired "$fwd_fastq_file" "$rev_fastq_file"
done

