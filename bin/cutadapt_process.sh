#!/usr/bin/env bash

set -euo pipefail

# Function to print help
usage() {
    echo "Usage: $0 -1 forward_read -2 reverse_read -m cutadapt_minlen -f fwd_primers_file -r rev_primers_file -s sequencer -e allowed_errors -c cores -o output"
    exit 1
}

# default settings
cores=1
allowed_errors=0 # allow no mismatches in the adapter sequence

# Parse command-line options
while getopts "1:2:m:f:r:s:e:c:o:" OPTION; do
    case $OPTION in
        1) forward_read=$OPTARG ;;
        2) reverse_read=$OPTARG ;;
        m) cutadapt_minlen=$OPTARG ;;
        f) fwd_primers_file=$OPTARG ;;
        r) rev_primers_file=$OPTARG ;;
        e) allowed_errors=$OPTARG ;;
        s) sequencer=$OPTARG ;;
        c) cores=$OPTARG ;;
        o) trimmed_demuxed_fastqs=$OPTARG ;;
        h) usage ;;
        ?) usage ;;
    esac
done

# Validate inputs
if [[ -z "$forward_read" || -z "$reverse_read" || -z "$cutadapt_minlen" || -z "$fwd_primers_file" || -z "$rev_primers_file" || -z "$sequencer" || -z "$trimmed_demuxed_fastqs" ]]; then
    usage
fi

# Check if input files exist
if [[ ! -f "$forward_read" ]]; then
    echo "Error: Forward read file ($forward_read) does not exist." >&2
    exit 2
fi

if [[ ! -f "$reverse_read" ]]; then
    echo "Error: Reverse read file ($reverse_read) does not exist." >&2
    exit 2
fi

# Intermediate directories for different stages
adapter_dimers="adapter_dimers"
no_adapter_dimers="no_adapter_dimers"
too_short="too_short"
fastp_output_preprocess="fastp_reports/preprocess"
fastp_output_postprocess="fastp_reports/postprocess"
trimmed_demuxed_unknown_fastqs="trimmed_demuxed_unknown_fastqs"

# Create directories
for dir in $adapter_dimers $no_adapter_dimers $too_short $fastp_output_preprocess $fastp_output_postprocess $trimmed_demuxed_unknown_fastqs $trimmed_demuxed_fastqs; do
    test -d $dir || mkdir -p $dir
done

cutadapt_json="cutadapt.json"

# Sample ID extraction
sample_id=$(basename "$forward_read" | sed 's/_R[12].*//')

# Step 0: Quality control before further processing
fastp \
    -i ${forward_read} \
    -I ${reverse_read} \
    -o ${adapter_dimers}/${sample_id}_R1.fastq.gz \
    -O ${adapter_dimers}/${sample_id}_R2.fastq.gz \
    --json ${fastp_output_preprocess}/${sample_id}_before_trim.json \
    --html ${fastp_output_preprocess}/${sample_id}_before_trim.html \
    --report_title "Before Trimming - ${sample_id}" \
    --detect_adapter_for_pe \
    --overrepresentation_analysis \
    --thread ${cores} \
    -Q \
    -L \
    -G \
    -A

# Step 1: Adapter trimming with relaxed criteria
cutadapt \
    --action=trim \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -e ${allowed_errors} \
    --no-indels \
    --minimum-length ${cutadapt_minlen} \
    -o ${adapter_dimers}/${sample_id}_dimer_R1.fastq.gz \
    -p ${adapter_dimers}/${sample_id}_dimer_R2.fastq.gz \
    --cores ${cores} \
    --untrimmed-output ${no_adapter_dimers}/${sample_id}_filtered_R1.fastq.gz \
    --untrimmed-paired-output ${no_adapter_dimers}/${sample_id}_filtered_R2.fastq.gz \
    --json=${cutadapt_json} \
    --compression-level=1 \
    --quiet \
    ${adapter_dimers}/${sample_id}_R1.fastq.gz \
    ${adapter_dimers}/${sample_id}_R2.fastq.gz > /dev/null

# Step 2: Quality control after adapter trimming
fastp \
    -i ${no_adapter_dimers}/${sample_id}_filtered_R1.fastq.gz \
    -I ${no_adapter_dimers}/${sample_id}_filtered_R2.fastq.gz \
    --json ${fastp_output_postprocess}/${sample_id}_after_trim.json \
    --html ${fastp_output_postprocess}/${sample_id}_after_trim.html \
    --report_title "After Trimming - ${sample_id}" \
    --detect_adapter_for_pe \
    --overrepresentation_analysis \
    --thread ${cores} \
    -Q \
    -L \
    -G \
    -A

# Step 3: Primer removal and demultiplexing
cutadapt \
    --action=trim \
    -g file:${fwd_primers_file} \
    -G file:${rev_primers_file} \
    --pair-adapters \
    -e ${allowed_errors} \
    --no-indels \
    --minimum-length ${cutadapt_minlen} \
    -o ${trimmed_demuxed_fastqs}/{name}_${sample_id}_trimmed_R1.fastq.gz \
    -p ${trimmed_demuxed_fastqs}/{name}_${sample_id}_trimmed_R2.fastq.gz \
    --untrimmed-output ${trimmed_demuxed_unknown_fastqs}/${sample_id}_unknown_R1.fastq.gz \
    --untrimmed-paired-output ${trimmed_demuxed_unknown_fastqs}/${sample_id}_unknown_R2.fastq.gz \
    --json=${cutadapt_json} \
    --compression-level=1 \
    --too-short-output ${too_short}/${sample_id}_too_short_R1.fastq.gz \
    --too-short-paired-output ${too_short}/${sample_id}_too_short_R2.fastq.gz \
    --quiet \
    ${no_adapter_dimers}/${sample_id}_filtered_R1.fastq.gz \
    ${no_adapter_dimers}/${sample_id}_filtered_R2.fastq.gz > /dev/null

# Step 4: Quality control on demultiplexed FASTQs
for fastq_file in ${trimmed_demuxed_fastqs}/*_trimmed_R1.fastq.gz; do
    base_name=$(basename $fastq_file _trimmed_R1.fastq.gz)
    fastp \
        -i ${trimmed_demuxed_fastqs}/${base_name}_trimmed_R1.fastq.gz \
        -I ${trimmed_demuxed_fastqs}/${base_name}_trimmed_R2.fastq.gz \
        --json ${fastp_output_postprocess}/${base_name}_demuxed.json \
        --html ${fastp_output_postprocess}/${base_name}_demuxed.html \
        --report_title "Demultiplexed - ${base_name}" \
        --detect_adapter_for_pe \
        --overrepresentation_analysis \
        --thread ${cores} \
        -Q \
        -L \
        -G \
        -A \
        -P 1 
done

# Extract read counts from cutadapt JSON
total_pairs=$(jq '.read_counts.input' ${cutadapt_json})
no_dimers=$(jq '.read_counts.filtered.discard_untrimmed' ${cutadapt_json})
too_short_count=$(jq '.read_counts.filtered.too_short' ${cutadapt_json})

# Count reads in each forward (F) file in trimmed_demuxed_unknown_fastqs
unknown_count=0
for file in ${trimmed_demuxed_unknown_fastqs}/*_unknown_R1.fastq.gz; do
    if [[ -f $file ]]; then
        read_count=$(zgrep -c ^@ $file) || zgrep_exit_status=$?
        zgrep_exit_status=${zgrep_exit_status:-0}
        if [[ $zgrep_exit_status -eq 1 ]]; then
            read_count=0
        elif [[ $zgrep_exit_status -ne 0 ]]; then
            echo "Error: could not count reads for $file"
            exit $zgrep_exit_status
        fi
        unknown_count=$((unknown_count + read_count))
    fi

    # unset the status variable
    unset zgrep_exit_status
done

# Summary of reads at different stages
printf "%s\t%s\n" "Input" ${total_pairs} > ${sample_id}.SAMPLEsummary.txt
printf "%s\t%s\n" "No Dimers" ${no_dimers} >> ${sample_id}.SAMPLEsummary.txt
printf "%s\t%s\n" "Too Short Reads" ${too_short_count} >> ${sample_id}.SAMPLEsummary.txt
printf "%s\t%s\n" "Unknown Reads" ${unknown_count} >> ${sample_id}.SAMPLEsummary.txt

# Count reads in each demultiplexed fastq file using zgrep
amplicon_counts=${sample_id}.AMPLICONsummary.txt
total_count=0

touch $amplicon_counts

while read -r amplicon_name; do
    fastq_file="${trimmed_demuxed_fastqs}/${amplicon_name}_${sample_id}_trimmed_R1.fastq.gz"
    if [[ -f $fastq_file ]]; then
        read_count=$(zgrep -c ^@ $fastq_file) || zgrep_exit_status=$?
        zgrep_exit_status=${zgrep_exit_status:-0}
        if [[ $zgrep_exit_status -eq 1 ]]; then
            read_count=0
        elif [[ $zgrep_exit_status -ne 0 ]]; then
            echo "Error: could not calculate amplicons for $fastq_file"
            exit $zgrep_exit_status
        fi
    else
        read_count=0
    fi

    total_count=$((total_count + read_count))
    echo -e "${amplicon_name}\t${read_count}" >> $amplicon_counts
    
    # unset the status variable
    unset zgrep_exit_status

done < <(grep '^>' $fwd_primers_file | sed 's/^>//')

printf "%s\t%s\n" "Amplicons" ${total_count} >> ${sample_id}.SAMPLEsummary.txt
