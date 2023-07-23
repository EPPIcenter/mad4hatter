#!/usr/bin/env bash
set -e

# Usage message
usage() {
  echo "Usage: $0 -a amplicon_table -f forward_primers -r reverse_primers
This script parses the mad4hatter amplicon table provided by the user and creates forward and reverse primer files from it.

Options:
    -h    Display this help and exit

  IN:
    -a    Amplicon table

  OUT:
    -f    Forward primers file
    -r    Reverse primers file" 1>&2
  exit 1
}

# Initialize parameters
fwd_primers_file=""
rev_primers_file=""
amplicon_info=""

# Parse options
while getopts "f:r:h:a:" o; do
  case "${o}" in
    a)
      amplicon_info=${OPTARG}
      ;;
    f)
      fwd_primers_file=${OPTARG}
      ;;
    r)
      rev_primers_file=${OPTARG}
      ;;
    h)
      usage
      ;;
    *)  # default case
      usage
      ;;
  esac
done
shift $((OPTIND-1))

# Check if both file paths were provided
if [ -z "${amplicon_info}" ] || [ -z "${fwd_primers_file}" ] || [ -z "${rev_primers_file}" ]; then
  usage
fi

# Create a temporary file to hold the forward and reverse adapters
adapters=$(mktemp)

# Delete the file on exit
trap "rm -f '$tmpfile'" EXIT

# Parse the amplicon_info file and create adapters.txt
cat $amplicon_info | awk 'NR==1 {
  for (i = 1; i <= NF; i++) {
    if ( $i == "amplicon" ) {amplicon=i}
    if ( $i == "fwd_primer" ) {fwd_primer=i}
    if ( $i == "rev_primer" ) {rev_primer=i}
  }
} NR>=1 {
  print $amplicon,$fwd_primer,$rev_primer
}' > ${adapters}

# Check if adapters.txt has 3 fields.
if [[ $(head -n 1 ${adapters} | awk '{print NF}') -ne 3 ]]; then 
  echo "ERROR: Must have 'fwd_primer' and 'rev_primer' in ${amplicon_info}!!!"
  exit 1 
fi

# Generate forward and reverse primer files
cut -d ' ' -f 1,2 $adapters | awk 'NR>1 { printf(">%s\n^%s\n", $1, $2) }' > $fwd_primers_file
cut -d ' ' -f 1,3 $adapters | awk 'NR>1 { printf(">%s\n^%s\n", $1, $2) }' > $rev_primers_file