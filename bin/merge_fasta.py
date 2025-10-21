#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import sys
import os


def parse_args_merge_fasta():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference_paths', type=str, required=True, nargs='+',
                        help='Paths to amplicon info tables for pools')
    parser.add_argument('--reference_output_path',
                        type=str, default='reference.fasta')
    return parser.parse_args()


def merge_fasta():
    args = parse_args_merge_fasta()
    # A dictionary to store unique sequences
    sequences = {}

    for file in args.reference_paths:
        if not os.path.isfile(file):
            print(f"File {file} does not exist. Skipping.")
            continue

        # Read sequences from the file
        for record in SeqIO.parse(file, "fasta"):
            seq_str = str(record.seq)
            if seq_str not in sequences:
                sequences[seq_str] = record.id
            else:
                print(f"Duplicate found: {record.id} (Skipping)")
    # Write the unique sequences to the output file
    with open(args.reference_output_path, "w") as output_handle:
        for seq, id in sequences.items():
            output_handle.write(f">{id}\n{seq}\n")


if __name__ == "__main__":
    merge_fasta()
