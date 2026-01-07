#!/usr/bin/env python3
import argparse
import os
import sys
import warnings
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter


def parse_args_merge_fasta():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--reference_paths', 
        type=str, 
        required=True, 
        nargs='+',
        help='Paths to amplicon info tables for pools'
        )
    parser.add_argument(
        '--reference_output_path',
        type=str, 
        default='reference.fasta'
        )
    parser.add_argument(
        '--wrap',
        type=int,
        default=80,
        help='FASTA line width (default: 80)'
    )
    return parser.parse_args()


def merge_fasta():
    args = parse_args_merge_fasta()
    # A dictionary to store unique sequences
    sequences = {}

    for file in args.reference_paths:
        if not os.path.isfile(file):
            warnings.warn(
                f"File {file} does not exist. Skipping.",
                category=UserWarning
            )
            continue

        # Read sequences from the file
        for record in SeqIO.parse(file, "fasta"):
            seq_str = str(record.seq)
            if seq_str not in sequences:
                sequences[seq_str] = record.id
            else:
                # print(f"Duplicate found: {record.id} (Skipping)", file=sys.stderr)
                warnings.warn(
                    f"Duplicate sequence found: {record.id} (Skipping)",
                    category=UserWarning
                )
    # Write the unique sequences to the output file
    if not sequences:
        sys.exit("ERROR: No valid FASTA records were found.")

    # Convert to SeqRecord objects
    records = sorted(
        (
            SeqRecord(
                Seq(seq),
                id=record_id,
                description=""
            )
            for seq, record_id in sequences.items()
        ),
        key=lambda r: r.id
    )

    # Write FASTA with 30 characters per line
    with open(args.reference_output_path, "w") as output_handle:
        writer = FastaWriter(output_handle, wrap=args.wrap)
        writer.write_file(records)


if __name__ == "__main__":
    merge_fasta()
