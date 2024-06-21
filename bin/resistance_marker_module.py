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


def extract_reference_per_marker(ref_marker_table, ref_sequences):
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
        return {
            'translatedCodonStart': codon_start,
            'CodonEnd': codon_end,
            'RefCodon': ref_codon
        }
    new_columns = ref_marker_table.apply(process_row, axis=1)
    new_columns_df = pd.DataFrame(
        new_columns.tolist(), index=ref_marker_table.index)

    ref_marker_table = pd.concat([ref_marker_table, new_columns_df], axis=1)

    # Translate RefCodon to amino acis
    ref_marker_table['RefAA'] = ref_marker_table['RefCodon'].apply(
        lambda x: translate(x))
    return ref_marker_table


def extract_mask_coordinates(mask_string):
    start, length = map(int, mask_string.strip('N').split('+'))
    start = start - 1
    end = start + length
    return list(range(start, end))


def get_locus_mask_coordinates(pseudo_cigar):
    # Get masking coordinates
    masking_pattern = r'\d+\+\d+N'
    masks = re.findall(masking_pattern, pseudo_cigar)
    mask_coordinates = []
    for mask in masks:
        mask_coordinates.extend(extract_mask_coordinates(mask))
    return set(mask_coordinates)


def create_masked_coordinate_table(allele_data):
    # Only process the first cigar from each locus
    cigar_per_locus = allele_data.groupby('Locus', as_index=False).first()

    cigar_per_locus['masked_coords'] = cigar_per_locus['PseudoCIGAR'].apply(
        get_locus_mask_coordinates)

    return cigar_per_locus[['Locus', 'masked_coords']]


def set_codon_with_indel_as_undetermined(codon, ref_codon):
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
        codon = set_codon_with_indel_as_undetermined(codon, aligned_ref_codon)
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


def parse_pseudo_cigar(string):
    """
    Parse a pseudo CIGAR string into a list of tuples (position, operation, value)
    The pseudo CIGAR string is a string representation of the differences between a reference sequence and a query sequence. The string
    is composed of a series of operations and values. The operations are "I" (insertion), "D" (deletion), and "N" (mask). Substitutions
    are also represented in the string as '[position][snp]'. For example, you could see "3C", where "3" is the position of the base
    and "C" is the base in the query sequence.

    :Example:
    '2N+45I=ACT7G' # Mask beginning at position 2 and extends 4 bases, 'ACT' insertion at position 5, 'G' substitution at position 7
    :param string: Pseudo CIGAR string
    :param orientation: "+" or "-"
    :param refseq_len: Length of the reference sequence
    :return: List of tuples (position, operation, value)
    """

    logging.debug(f"Entering `parse_pseudo_cigar` with string={
                  string}")

    if not isinstance(string, str):
        logging.error(f"Expected a string for `parse_pseudo_cigar` but got {
                      type(string)}. Value: {string}")
        raise TypeError(f"Expected a string for `parse_pseudo_cigar` but got {
                        type(string)}. Value: {string}")

    string = string.replace("=", "")  # Remove the "=" separator
    logging.debug(f"Removed '=' separator. Modified string={string}")

    tuples = []

    pattern = re.compile(r'(\d+)([ID\+]?)\+?(\d+|[ACGTN]+)?')
    matches = pattern.findall(string)
    logging.debug(f"Pattern found {len(matches)} matches in the string")

    for match in matches:
        position, operation, value = match
        position = int(position) - 1  # use base 0 indexing
        logging.debug(f"Processing match: original position={match[0]}, adjusted position={
                      position}, operation={operation}, value={value}")

        # Handle mask
        if operation == "+":
            logging.debug("Operation is a mask. No changes applied.")
            pass  # do nothing for now
        else:
            tuples.append((position, None if operation ==
                          '' else operation, value))
            logging.debug(f"Added tuple to result: {
                          (position, None if operation == '' else operation, value)}")

    logging.debug(f"Exiting `parse_pseudo_cigar` with {
                  len(tuples)} tuples generated")
    return tuples


def extract_mutations_from_df(df, ref_sequences):
    results = []
    transtab = str.maketrans("TACG", "ATGC")
    for index, row in df.iterrows():
        pseudo_cigar = row['PseudoCIGAR']
        locus = row['Locus']
        strand = row['strand']

        if pseudo_cigar == ".":
            continue

        # returns 0-based position
        changes = parse_pseudo_cigar(pseudo_cigar)

        ref_seq = ref_sequences[locus].seq
        for pos, op, alt in changes:
            ref = ref_seq[int(pos)]
            if strand == '-':
                alt = alt.translate(transtab)
                ref = ref.translate(transtab)
                pos = len(ref_seq) - pos - 1

            results.append({
                'PseudoCIGAR': pseudo_cigar,
                'Position': pos,
                'Alt': alt,
                'Ref': ref,
            })

    return pd.DataFrame(results)


def generate_resmarker_table(asv_table, allele_data):
    # Get resmarker codon
    resmarker_per_unique_asv = asv_table.apply(
        apply_process_resmarker, axis=1)
    # Merge the results back to the original DataFrame
    resmarker_per_unique_asv = pd.concat(
        [asv_table, resmarker_per_unique_asv], axis=1)

    # Add on additional columns
    resmarker_per_unique_asv['AA'] = resmarker_per_unique_asv['Codon'].apply(
        lambda x: translate(x))
    resmarker_per_unique_asv['CodonRefAlt'] = np.where(
        resmarker_per_unique_asv['RefCodon'] == resmarker_per_unique_asv['Codon'], 'REF', 'ALT')
    resmarker_per_unique_asv['AARefAlt'] = np.where(
        resmarker_per_unique_asv['RefAA'] == resmarker_per_unique_asv['AA'], 'REF', 'ALT')

    # Merge back with samples
    resmarker_data = resmarker_per_unique_asv.merge(
        allele_data[['SampleID', 'Locus', 'ASV', 'PseudoCIGAR', 'Reads', 'CodonStart']], on=['Locus', 'ASV', 'PseudoCIGAR', 'CodonStart'])

    # Condense down repeat codons
    resmarker_data = resmarker_data.groupby(['SampleID', 'Locus', 'RefCodon', 'Codon', 'gene_id',
                                             'translatedCodonStart', 'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt', 'strand']).Reads.sum().reset_index()

    # Ensure reads is integer
    resmarker_data['Reads'] = resmarker_data['Reads'].astype(int)
    return resmarker_data


def main(args):

    logging.debug(f"------ Start of `main` ------")
    # READ DATA
    allele_data = pd.read_csv(args.allele_data_path, sep='\t')
    aligned_asv_data = pd.read_csv(args.aligned_asv_table_path, sep='\t')
    res_markers_info = pd.read_csv(args.res_markers_info_path, sep='\t')
    ref_sequences = SeqIO.to_dict(SeqIO.parse(args.refseq_path, 'fasta'))

    # GET MASKED COORDS - these are 1-based
    masked_coords_per_locus = create_masked_coordinate_table(allele_data)

    # GET REFERENCE FOR MARKERS
    ref_marker_table = extract_reference_per_marker(
        res_markers_info[['gene_id', 'strand', 'Locus', 'CodonStart']], ref_sequences)

    # MERGE DATA
    # Merge allele data with resmarkers, only keeping the loci we need and markers with calls
    allele_data_per_resmarker = allele_data.merge(
        res_markers_info, on='Locus', how='inner')
    # Merge with alignment data
    allele_data_per_resmarker = allele_data_per_resmarker.merge(
        aligned_asv_data, left_on=['SampleID', 'Locus', 'ASV'],
        right_on=['sampleID', 'refid', 'asv'], how='left')

    # GET TABLE TO EXTRACT INFO FROM
    # Compile info for unique asv and markers
    required_columns = ['Locus', 'ASV', 'PseudoCIGAR',
                        'strand', 'gene_id', 'CodonStart', 'hapseq', 'refseq']
    unique_asvs_to_process = allele_data_per_resmarker[required_columns].drop_duplicates(
    )
    unique_asvs_to_process = unique_asvs_to_process.merge(ref_marker_table, on=['gene_id', 'strand',
                                                                                'Locus', 'CodonStart'], how='left')
    unique_asvs_to_process = unique_asvs_to_process.merge(
        masked_coords_per_locus, on='Locus', how='left')

    # Get resmarker table
    resmarker_data = generate_resmarker_table(
        unique_asvs_to_process, allele_data_per_resmarker)

    # Create microhaplotype version
    mhap_table = generate_mhap_table(resmarker_data)

    # Get all mutations
    mutations_df = extract_mutations_from_df(
        unique_asvs_to_process, ref_sequences)
    all_mutations = allele_data.merge(mutations_df, on='PseudoCIGAR')[
        ['SampleID', 'Locus', 'PseudoCIGAR', 'Position', 'Alt', 'Ref', 'Reads']]
    all_mutations = all_mutations.groupby(['SampleID', 'Locus', 'PseudoCIGAR',
                                           'Position', 'Alt', 'Ref']).Reads.sum().reset_index()

    resmarker_data.to_csv('resmarker_table.txt', sep='\t', index=False)
    mhap_table.to_csv('resmarker_microhap_table.txt', sep='\t', index=False)
    all_mutations.to_csv('discovery_table.txt', sep='\t', index=False)


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

# TODO: Add logging back in with format
# TODO: change AARefAlt to be ref alt instead of true false
# TODO: format grouping of the outputs
# TODO: turn these into classes
# TODO: decide on position to show : either strand and position realtive to start or just the one in the other table (or both)
# TODO: Decide on whether to have geneID split up
# TODO: Check if masking works 0 based?
print(datetime.now() - startTime)
