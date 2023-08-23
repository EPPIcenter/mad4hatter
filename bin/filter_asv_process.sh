#!/usr/bin/env bash

# Function to print help
usage() {
    echo "Usage: $0 -i alignments.txt -o alignments.filtered.txt -t filter_threshold"
    exit 1
}

# Set default values
colname="score"

# Parse command-line options
while getopts "i:o:t:" OPTION
do
    case $OPTION in
        i)
            infile=$OPTARG
            ;;
        o)
            outfile=$OPTARG
            ;;
        t)
            threshold=$OPTARG
            ;;
        ?)
            echo "Caught unknown: $OPTION"
            usage
            ;;
    esac
done

echo "infile: $infile"
echo "outfile: $outfile"
echo "threshold: $theshold"

# Validate inputs
if [[ -z "$infile" || -z "$outfile" || -z "$threshold" ]]; then
    usage
fi

# Find the column number for 'score'
col_num=$(head -1 "$infile" | tr '\t' '\n' | nl | grep -iw $colname | awk '{print $1}')

if [[ -z "$col_num" ]]; then
    echo "Error: 'score' column not found."
    exit 2
fi

# Use awk to filter based on the threshold
awk -v col="$col_num" -v thresh="$threshold" 'NR == 1 || $col > thresh' "$infile" > "$outfile"

echo "Filtered data saved to $outfile"