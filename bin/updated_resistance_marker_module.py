#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from Bio import SeqIO
from Bio.Seq import translate, Seq
import argparse
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import numpy as np
import json
import logging
from logging.handlers import RotatingFileHandler
import traceback
from datetime import datetime
startTime = datetime.now()


def create_reference_table(ref_table, ref_sequences):
    def process_row(row):
        locus_seq = ref_sequences[row['Locus']].seq
        if row['strand'] == '-':
            refseq_len = len(locus_seq)
            refseq_rc = locus_seq.reverse_complement()
            codon_start = refseq_len - \
                row['CodonStart'] - 2  # 0-based indexing
            codon_end = refseq_len - row['CodonStart'] + 1
            ref_codon = str(refseq_rc[codon_start: codon_end])
        else:
            codon_start = row['CodonStart'] - 1
            codon_end = row['CodonStart'] + 2
            ref_codon = str(locus_seq[codon_start: codon_end])
        return codon_start, codon_end, ref_codon

    # Apply the process_row function to each row in ref_table
    results = ref_table.apply(lambda row: process_row(row), axis=1)

    # Split the results into separate columns
    ref_table[['translatedCodonStart', 'CodonEnd', 'RefCodon']] = pd.DataFrame(
        results.tolist(), index=ref_table.index)

    # Translate RefCodon to RefAA using BioPython's translate
    ref_table['RefAA'] = ref_table['RefCodon'].apply(
        lambda x: translate(x))

    return ref_table


def get_mask_coordinates(pseudo_cigar):
    # Get masking coordinates
    masking_pattern = r'\d+\+\d+N'
    masks = re.findall(masking_pattern, pseudo_cigar)
    mask_coordinates = []
    for mask in masks:
        mask_coordinates.extend(extract_mask_coordinates(mask))
    return set(mask_coordinates)


def extract_mask_coordinates(mask):
    start, length = map(int, mask.strip('N').split('+'))
    end = start + length
    return list(range(start, end))


def create_masked_coordinate_table(allele_data):
    # TODO: double check the basing of this is correct
    # Group by 'locus' and get the first 'pseudo_cigar' for each group
    cigar_per_locus = allele_data.groupby('Locus', as_index=False).first()

    # Apply 'get_mask_coordinates' to each row's 'pseudo_cigar'
    cigar_per_locus['masked_coords'] = cigar_per_locus['PseudoCIGAR'].apply(
        get_mask_coordinates)

    return cigar_per_locus[['Locus', 'masked_coords']]


def set_codon_as_undetermined(codon, ref_codon):
    if '-' in codon or '-' in ref_codon:
        return 'X'
    else:
        return codon


def extract_codon_from_asv(reference, sequence, mask_coordinates, start_position):
    preceding_indel_present = False
    # No indels
    if '-' not in reference and '-' not in sequence:
        codon = sequence[start_position:start_position+3]
    # Indels
    else:
        index = 0
        while index < start_position:
            ref_base = reference[index]
            seq_base = sequence[index]
            if ref_base == '-':
                start_position += 1
                if index + 1 not in mask_coordinates:
                    preceding_indel_present = True
            if seq_base == '-':
                if index + 1 not in mask_coordinates:
                    preceding_indel_present = True
            index += 1
        codon = sequence[start_position:start_position+3]
        aligned_ref_codon = reference[start_position:start_position+3]
        codon = set_codon_as_undetermined(codon, aligned_ref_codon)
    return codon, preceding_indel_present


def process_resmarker(reference, sequence, start_position, pseudo_cigar, ref_codon, mask_coordinates, strand):
    if pseudo_cigar == '.':
        codon = ref_codon
        codon_masked = False
        preceding_indel_present = False
    else:
        if strand == '-':
            sequence = Seq(sequence).reverse_complement()
            reference = Seq(reference).reverse_complement()
        codon_coordinates = list(range(start_position, start_position+3))

        # Check if the codon is masked
        codon_masked = bool(mask_coordinates.intersection(codon_coordinates))

        # Get codon from asv
        codon, preceding_indel_present = extract_codon_from_asv(
            reference, sequence, mask_coordinates, start_position
        )

    return codon, codon_masked, preceding_indel_present


def apply_process_resmarker(row):
    codon, codon_masked, preceding_indel_present = process_resmarker(
        row['refseq'], row['hapseq'], row['translatedCodonStart'], row['PseudoCIGAR'], row['RefCodon'], row['masked_coords'], row['strand']
    )
    return pd.Series({'Codon': codon, 'codon_masked': codon_masked, 'preceding_indel_present': preceding_indel_present})


def generate_mhap_table(resmarker_table):
    resmarker_table = resmarker_table.sort_values(by='translatedCodonStart')
    mhap_table = resmarker_table.groupby(['SampleID', 'Locus', 'Reads']).agg(
        MicrohapIndex=('translatedCodonStart',
                       lambda x: '/'.join(map(str, x))),
        RefMicrohap=('RefAA', lambda x: '/'.join(map(str, x))),
        Microhaplotype=('AA', lambda x: '/'.join(map(str, x)))
    ).reset_index()
    mhap_table['MicrohapRefAlt'] = np.where(
        mhap_table['RefMicrohap'] == mhap_table['Microhaplotype'], 'REF', 'ALT')
    return mhap_table


def main(args):

    logging.debug(f"------ Start of `main` ------")
    # READ DATA
    allele_data = pd.read_csv(args.allele_data_path, sep='\t')
    aligned_asv_data = pd.read_csv(args.aligned_asv_table_path, sep='\t')
    res_markers_info = pd.read_csv(args.res_markers_info_path, sep='\t')
    ref_sequences = SeqIO.to_dict(SeqIO.parse(args.refseq_path, 'fasta'))

    # MERGE DATA
    # Merge allele data with resmarkers, only keeping the loci we need and markers with calls
    allele_data = allele_data.merge(
        res_markers_info, on='Locus', how='inner')
    # Merge with alignment data
    allele_data = allele_data.merge(
        aligned_asv_data, left_on=['SampleID', 'Locus', 'ASV'],
        right_on=['sampleID', 'refid', 'asv'], how='left')

    # GET MASKED COORDS
    masked_coords_per_locus = create_masked_coordinate_table(allele_data)

    # GET REFERENCE FOR MARKERS
    # TODO: could use the ref sequence from the alignments for this
    ref_table = allele_data[['gene_id', 'strand',
                             'Locus', 'CodonStart']].drop_duplicates()
    ref_table = create_reference_table(ref_table, ref_sequences)

    # GET TABLE TO EXTRACT INFO FROM
    # Compile info for unique asv and markers
    required_columns = ['Locus', 'ASV', 'PseudoCIGAR',
                        'strand', 'gene_id', 'CodonStart', 'hapseq', 'refseq']
    unique_asvs_to_process = allele_data[required_columns].drop_duplicates()
    # Add on the ref codon - could change this to left
    unique_asvs_to_process = unique_asvs_to_process.merge(ref_table, on=['gene_id', 'strand',
                                                                         'Locus', 'CodonStart'], how='outer')
    # Add on masked coordinates - could change this to left
    unique_asvs_to_process = unique_asvs_to_process.merge(
        masked_coords_per_locus, on='Locus', how='outer')
    # Get resmarker codon
    unique_asvs_to_process_results = unique_asvs_to_process.apply(
        apply_process_resmarker, axis=1)
    # Merge the results back to the original DataFrame
    unique_asvs_to_process = unique_asvs_to_process.reset_index(drop=True)
    unique_asvs_to_process = pd.concat(
        [unique_asvs_to_process, unique_asvs_to_process_results], axis=1)

    # Add on additional information
    unique_asvs_to_process['AA'] = unique_asvs_to_process['Codon'].apply(
        lambda x: translate(x))
    unique_asvs_to_process['CodonRefAlt'] = np.where(
        unique_asvs_to_process['RefCodon'] == unique_asvs_to_process['Codon'], 'REF', 'ALT')
    unique_asvs_to_process['AARefAlt'] = np.where(
        unique_asvs_to_process['RefAA'] == unique_asvs_to_process['AA'], 'REF', 'ALT')
    # TODO: Make GeneID  Gene    CodonID
    # Merge back with samples
    resmarker_data = unique_asvs_to_process.merge(
        allele_data[['SampleID', 'Locus', 'ASV', 'PseudoCIGAR', 'Reads', 'CodonStart']], on=['Locus', 'ASV', 'PseudoCIGAR', 'CodonStart'])
    # Condense down repeat codons
    resmarker_data = resmarker_data.groupby(['SampleID', 'Locus', 'RefCodon', 'Codon', 'gene_id',
                                             'translatedCodonStart', 'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt', 'strand']).Reads.sum().reset_index()

    # Ensure reads is integer
    resmarker_data['Reads'] = resmarker_data['Reads'].astype(int)
    resmarker_data.to_csv('full_merged_data.csv', sep='\t', index=False)

    # Create microhaplotype version
    mhap_table = generate_mhap_table(resmarker_data)

    mhap_table.to_csv('mhaps.txt')


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(
            description="Process some input files.")

        parser.add_argument("--allele_data_path",
                            required=True, help="Path to allele_data.txt")
        parser.add_argument("--aligned_asv_table_path",
                            required=True, help="Path to the aligned asv table")
        parser.add_argument("--res_markers_info_path",
                            required=True, help="Path to resistance marker table")
        parser.add_argument("--refseq_path", required=True,
                            help="Path to reference sequences [fasta]")
        parser.add_argument("--n-cores", type=int, default=1,
                            help="Number of cores to use")
        parser.add_argument("--log-level", default="INFO", choices=[
                            "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Set the logging level")
        parser.add_argument("--log-file", type=str,
                            help="Write logs to a file instead of the console")
        parser.add_argument("--log-max-size", type=int,
                            default=5, help="Maximum log file size in MB")
        parser.add_argument("--log-backups", type=int,
                            default=1, help="Number of log files to keep")

        args = parser.parse_args()

        numeric_level = getattr(logging, args.log_level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f"Invalid log level: {args.log_level}")

        if args.log_file:
            # Use rotating log files
            file_handler = RotatingFileHandler(
                args.log_file, maxBytes=args.log_max_size * 1024 * 1024, backupCount=args.log_backups)
            logging.basicConfig(
                level=numeric_level, format="%(asctime)s - %(levelname)s - %(message)s", handlers=[file_handler])
        else:
            logging.basicConfig(
                level=numeric_level, format="%(asctime)s - %(levelname)s - %(message)s")

        logging.info("Program start.")
        main(args)
    except Exception as e:
        # logging.error(f"An error of type {
        #   type(e).__name__} occurred. Message: {e}")
        logging.error(traceback.format_exc())  # Log the full traceback
        raise e
    finally:
        logging.info("Program complete.")
        logging.shutdown()

# Plan
# Load in allele data - DONE
# Create table of the masked coordinates for each locus - DONE
# Merge allele data with alignment data - DONE
# Merge the above table with the resmarker info - DONE
# Make a table of the reference for eachmarker - DONE
# For each unique asv and locus extract the corresponding resmarker - DONE
# Build the table

# TODO: Change so unique ASVs made first and then the other data merged on
# TODO: Add logging back in with format
# TODO: Add warning if there are lines in the resmarker file that are not covered by amplicon. make sure the codon is at least 3 bases from the end. Compare to info file
# TODO: make sure order groupby is correct
# TODO: Add on the new mutations
# TODO: change AARefAlt to be ref alt instead of true false
# TODO: turn these into classes
# TODO: groupby markers and codons and sum reads
print(datetime.now() - startTime)
