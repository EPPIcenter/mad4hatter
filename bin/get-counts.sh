#!/usr/bin/env bash

# Initialize default values
SPIKEIN_CSV=""
DEMUX_DIR=""
OUTPUT_DIR=""
THREADS=4

# Parse command-line options
while getopts "c:d:o:t:" opt; do
  case $opt in
    c) SPIKEIN_CSV="$OPTARG";;
    d) DEMUX_DIR="$OPTARG";;
    o) OUTPUT_DIR="$OPTARG";;
    t) THREADS="$OPTARG";;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1;;
  esac
done

test -z $SPIKEIN_CSV && { echo "Spikein csv is missing!!!"; exit 1; }

# Get sample name
SAMPLE=$(basename "$DEMUX_DIR")

# Define other derived variables
SPIKEIN_FASTA="${OUTPUT_DIR}/spikeins.fasta"
BWA_INDEX="${OUTPUT_DIR}/spikeins_index"
FINAL_COUNTS="${OUTPUT_DIR}/${SAMPLE}_final_spikein_counts.csv"
MULTI_MAP_COUNTS="${OUTPUT_DIR}/${SAMPLE}_multi_map_spikein_counts.csv"

# Ensure output directory exists
mkdir -p "${OUTPUT_DIR}"

# Header for the consolidated counts file
echo "SampleID,SpikeinID,Count" > "$FINAL_COUNTS"
echo "SampleID,SpikeinID,MultiMappedCount" > "$MULTI_MAP_COUNTS"

# Step 1: Convert the CSV file of spike-in sequences to FASTA format
awk -F, 'NR > 1 {print ">" $1 "\n" $2}' "$SPIKEIN_CSV" > "$SPIKEIN_FASTA"

# Step 2: Index the FASTA file with BWA
bwa index -p "$BWA_INDEX" "$SPIKEIN_FASTA"

ALIGNMENT_SAM="${OUTPUT_DIR}/${SAMPLE}_spikeins_aligned.sam" # outputs
SORTED_BAM="${OUTPUT_DIR}/${SAMPLE}_spikeins_aligned_sorted.bam" # outputs
FASTQ_R1="${DEMUX_DIR}/SDSI_${SAMPLE}_R1.fastq.gz" # inputs
FASTQ_R2="${DEMUX_DIR}/SDSI_${SAMPLE}_R2.fastq.gz" # inputs
SAMPLE_ALIGNMENT_STATS="${OUTPUT_DIR}/${SAMPLE}_alignment_stats.txt"

if [[ -f "$FASTQ_R1" && -f "$FASTQ_R2" ]]; then
    bwa mem -t "$THREADS" "$BWA_INDEX" "$FASTQ_R1" "$FASTQ_R2" > "$ALIGNMENT_SAM"
    samtools view -@ "$THREADS" -bS "$ALIGNMENT_SAM" | samtools sort -o "$SORTED_BAM"
    samtools index -@ "$THREADS" "$SORTED_BAM"
    samtools idxstats -@ "$THREADS" "$SORTED_BAM" > "${SAMPLE_ALIGNMENT_STATS}"
    awk -v sample="$SAMPLE" 'BEGIN{OFS=","} {print sample, $1, $3}' "${SAMPLE_ALIGNMENT_STATS}" >> "$FINAL_COUNTS"
    samtools view -@ "$THREADS" -q 30 -F 4 "$SORTED_BAM" | awk -v sample="$SAMPLE" '{print sample, $3}' | sort | uniq -c | awk '{print $2, $3, $1}' >> "$MULTI_MAP_COUNTS"

else
    echo "SDSI FASTQ files for $SAMPLE not found."
    exit 1
fi

echo "Consolidated counts saved to $FINAL_COUNTS"
