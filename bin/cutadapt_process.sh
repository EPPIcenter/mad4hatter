#!/usr/bin/env bash

set -e

# Function to print help
usage() {
    echo "Usage: $0 -1 forward_read -2 reverse_read -m cutadapt_minlen -f fwd_primers_file -r rev_primers_file -s sequencer -e allowed_errors -c cores -o output"
    exit 1
}

# default settings
cores=1
allowed_errors=0 # allow no mismatches in the adapter sequence

# Parse command-line options
while getopts "1:2:m:f:r:q:h:s:e:c:o:" OPTION
do
    case $OPTION in
        1)
            forward_read=$OPTARG
            ;;
        2)
            reverse_read=$OPTARG
            ;;
        m)
            cutadapt_minlen=$OPTARG
            ;;
        f)
            fwd_primers_file=$OPTARG
            ;;
        r)
            rev_primers_file=$OPTARG
            ;;
        e)
			allowed_errors=$OPTARG
			;;
        s)
            sequencer=$OPTARG
            ;;
        c)
			cores=$OPTARG
			;;
		o)
			trimmed_demuxed_fastqs=$OPTARG
			;;
        h)
            usage
            ;;
        ?)
            usage
            ;;
    esac
done

# Validate inputs
if [[ -z "$forward_read" || -z "$reverse_read" || -z "$cutadapt_minlen" || -z "$fwd_primers_file" || -z "$rev_primers_file" || -z "$sequencer" || -z "$trimmed_demuxed_fastqs" ]]; then
    usage
fi

# this is the directory that will contain forward and reverse
# reads that have been quality filtered and their primers removed
test -d $trimmed_demuxed_fastqs || mkdir -p $trimmed_demuxed_fastqs

# The remaining reads that could not be demultiplexed
# trimmed_demuxed_unknown_fastqs=$(mktemp -d)
trimmed_demuxed_unknown_fastqs="trimmed_demuxed_unknown_fastqs"
test -d $trimmed_demuxed_unknown_fastqs || mkdir -p $trimmed_demuxed_unknown_fastqs

# Create intermediate directories and files
# adapter_dimers=$(mktemp -d)
adapter_dimers="adapter_dimers"
test -d $adapter_dimers || mkdir -p $adapter_dimers

# no_adapter_dimers=$(mktemp -d)
no_adapter_dimers="no_adapter_dimers"
test -d $no_adapter_dimers || mkdir -p $no_adapter_dimers

# cutadapt_json=$(mktemp)
cutadapt_json="cutadapt.json"

# Delete the directories on exit if specified on the command line
trap "rm -rf '$adapter_dimers' '$no_adapter_dimers'" EXIT

# get the sample id
sample_id=$(basename "$forward_read" | sed 's/_R[12].*//')

# Remove all adapter dimers
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
    ${forward_read} \
    ${reverse_read} > /dev/null


# Get the total number of reads with adapters
total_pairs=$(jq '.read_counts.input' ${cutadapt_json})
no_dimers=$(jq '.read_counts.filtered.discard_untrimmed' ${cutadapt_json})
printf "%s\t%s\n" "Input" ${total_pairs} > ${sample_id}.SAMPLEsummary.txt
printf "%s\t%s\n" "No Dimers" ${no_dimers} >> ${sample_id}.SAMPLEsummary.txt

if [ "$sequencer" == "miseq" ]; then
    qualfilter="--trim-n -q 10"
else
    qualfilter="--nextseq-trim=20"
fi

cutadapt \
    --action=trim \
    -g file:${fwd_primers_file} \
    -G file:${rev_primers_file} \
    --pair-adapters \
    -e ${allowed_errors} \
    --no-indels \
    ${qualfilter} \
    --minimum-length ${cutadapt_minlen} \
    -o ${trimmed_demuxed_fastqs}/{name}_${sample_id}_trimmed_R1.fastq.gz \
    -p ${trimmed_demuxed_fastqs}/{name}_${sample_id}_trimmed_R2.fastq.gz \
    --untrimmed-output ${trimmed_demuxed_unknown_fastqs}/${sample_id}_unknown_R1.fastq.gz \
    --untrimmed-paired-output ${trimmed_demuxed_unknown_fastqs}/${sample_id}_unknown_R2.fastq.gz \
    --json=${cutadapt_json} \
    --compression-level=1 \
    --quiet \
    ${no_adapter_dimers}/${sample_id}_filtered_R1.fastq.gz \
    ${no_adapter_dimers}/${sample_id}_filtered_R2.fastq.gz > /dev/null


# Count reads in each demultiplexed fastq file using zgrep
amplicon_counts=${sample_id}.AMPLICONsummary.txt
sample_counts=${sample_id}.SAMPLEsummary.txt
total_count=0

touch $amplicon_counts $sample_counts

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

printf "%s\t%s\n" "Amplicons" ${total_count} >> ${sample_counts}