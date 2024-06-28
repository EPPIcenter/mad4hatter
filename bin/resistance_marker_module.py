#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from Bio import SeqIO
from Bio.Seq import translate, Seq
import argparse
import re
from functools import partial
import numpy as np
import json
import logging
from logging.handlers import RotatingFileHandler
import traceback
from datetime import datetime
startTime = datetime.now()


class ProcessMarkerInfo:
    def __init__(self, ref_sequences):
        self.ref_sequences = ref_sequences

    def extract_info_for_markers(self, res_marker_table):
        marker_info = res_marker_table.apply(
            self.extract_marker_info, axis=1)

        res_marker_table = pd.DataFrame(marker_info.tolist())

        # Translate RefCodon to amino acid
        res_marker_table['RefAA'] = res_marker_table['RefCodon'].apply(
            lambda x: translate(x))
        return res_marker_table

    def extract_marker_info(self, row):
        locus = row['Locus']
        strand = row['strand']
        codon_start = row['CodonStart']
        ref_seq = self.ref_sequences[locus].seq

        if strand == '-':
            refseq_len = len(ref_seq)
            ref_seq = ref_seq.reverse_complement()
            relative_codon_start = refseq_len - codon_start - 2
        else:
            relative_codon_start = codon_start - 1

        codon_end = relative_codon_start + 3
        ref_codon = str(ref_seq[relative_codon_start: codon_end])
        return {
            'GeneID': row.GeneID,
            'Gene': row.Gene,
            'CodonID': row.CodonID,
            'strand': strand,
            'Locus': locus,
            'relativeCodonStart': relative_codon_start,
            'CodonStart': codon_start,
            'CodonEnd': codon_end,
            'RefCodon': ref_codon
        }


class MaskCoordinatesExtractor:
    @staticmethod
    def mask_str_to_coordinates(mask_string):
        start, length = map(int, mask_string.strip('N').split('+'))
        start = start - 1
        end = start + length
        return list(range(start, end))

    @staticmethod
    def cigar_to_mask_coordinates(pseudo_cigar):
        # Get masking coordinates
        masking_pattern = r'\d+\+\d+N'
        masks = re.findall(masking_pattern, pseudo_cigar)
        mask_coordinates = []
        for mask in masks:
            mask_coordinates.extend(
                MaskCoordinatesExtractor.mask_str_to_coordinates(mask))
        return set(mask_coordinates)

    @staticmethod
    def create_masked_coordinate_table(allele_data):
        # Only process the first cigar from each locus
        cigar_per_locus = allele_data.groupby('Locus', as_index=False).first()

        cigar_per_locus['masked_coords'] = cigar_per_locus['PseudoCIGAR'].apply(
            MaskCoordinatesExtractor.cigar_to_mask_coordinates)

        return cigar_per_locus[['Locus', 'masked_coords']]


class CodonProcessor:
    @staticmethod
    def set_codon_with_indel_as_undetermined(codon, ref_codon):
        if '-' in codon or '-' in ref_codon:
            return 'X'
        else:
            return codon

    @staticmethod
    def extract_codon_from_asv(reference, asv, mask_coordinates, start_position, end_position):
        preceding_indel_present = False
        # No indels
        if '-' not in reference and '-' not in asv:
            codon = asv[start_position:end_position]
        # Indels
        else:
            index = 0
            while index < start_position:
                ref_base = reference[index]
                seq_base = asv[index]
                if ref_base == '-':
                    start_position += 1
                    if index + 1 not in mask_coordinates:
                        preceding_indel_present = True
                if seq_base == '-':
                    if index + 1 not in mask_coordinates:
                        preceding_indel_present = True
                index += 1
            codon = asv[start_position:end_position]
            aligned_ref_codon = reference[start_position:end_position]
            codon = CodonProcessor.set_codon_with_indel_as_undetermined(
                codon, aligned_ref_codon)
        return codon, preceding_indel_present

    @staticmethod
    def get_codon_info_for_marker(marker_info):
        reference = marker_info['refseq']
        sequence = marker_info['hapseq']
        start_position = marker_info['relativeCodonStart']
        pseudo_cigar = marker_info['PseudoCIGAR']
        ref_codon = marker_info['RefCodon']
        mask_coordinates = marker_info['masked_coords']
        strand = marker_info['strand']
        end_position = marker_info['CodonEnd']
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
            codon_masked = bool(
                mask_coordinates.intersection(codon_coordinates))

            # Get codon from asv
            codon, preceding_indel_present = CodonProcessor.extract_codon_from_asv(
                reference, sequence, mask_coordinates, start_position, end_position
            )
        marker_with_codon_info = (pd.concat([marker_info, pd.Series(
            {'Codon': codon, 'codon_masked': codon_masked, 'preceding_indel_present': preceding_indel_present})]))
        return marker_with_codon_info


def pseudo_cigar_to_mutations(string):
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
        changes = pseudo_cigar_to_mutations(pseudo_cigar)

        ref_seq = ref_sequences[locus].seq
        for pos, op, alt in changes:
            ref = ref_seq[int(pos)]
            if strand == '-':
                alt = alt.translate(transtab)
                ref = ref.translate(transtab)
                pos = len(ref_seq) - pos - 1

            results.append({
                'GeneID': row.GeneID,
                'Gene': row.Gene,
                'Locus': locus,
                'PseudoCIGAR': pseudo_cigar,
                'Position': pos,
                'Alt': alt,
                'Ref': ref,
            })

    return pd.DataFrame(results)


def generate_mhap_table(resmarker_table):
    resmarker_table = resmarker_table.sort_values(by='relativeCodonStart')
    mhap_table = resmarker_table.groupby(['SampleID', 'Locus', 'GeneID', 'Gene', 'PseudoCIGAR', 'Reads']).agg(
        MicrohapIndex=('CodonID',
                       lambda x: '/'.join(map(str, x))),
        RefMicrohap=('RefAA', lambda x: '/'.join(map(str, x))),
        Microhaplotype=('AA', lambda x: '/'.join(map(str, x)))
    ).reset_index()
    mhap_table['MicrohapRefAlt'] = np.where(
        mhap_table['RefMicrohap'] == mhap_table['Microhaplotype'], 'REF', 'ALT')
    return mhap_table


def generate_all_mutations_table(asv_table, allele_data, ref_sequences):
    # Get all mutations
    mutations_df = extract_mutations_from_df(
        asv_table, ref_sequences)
    all_mutations = allele_data.merge(mutations_df, on=['Locus', 'PseudoCIGAR'])[
        ['SampleID', 'Locus', 'GeneID', 'Gene', 'PseudoCIGAR', 'Position', 'Alt', 'Ref', 'Reads']]
    all_mutations = all_mutations.groupby(
        ['SampleID', 'Locus', 'GeneID', 'Gene', 'Position', 'Alt', 'Ref']).Reads.sum().reset_index()
    return all_mutations


def main(args):

    logging.debug(f"------ Start of `main` ------")

    # READ DATA
    allele_data = pd.read_csv(args.allele_data_path, sep='\t')
    aligned_asv_data = pd.read_csv(args.aligned_asv_table_path, sep='\t')
    res_markers_info = pd.read_csv(args.res_markers_info_path, dtype={
                                   'GeneID': str}, sep='\t')
    ref_sequences = SeqIO.to_dict(SeqIO.parse(args.refseq_path, 'fasta'))

    # MARKER INFO - Extract codon start relative to strand and reference information
    marker_info_extractor = ProcessMarkerInfo(ref_sequences)
    processed_res_markers_info = marker_info_extractor.extract_info_for_markers(
        res_markers_info)

    # GET MASKED COORDS
    mask_extractor = MaskCoordinatesExtractor()
    masked_coords_table = mask_extractor.create_masked_coordinate_table(
        allele_data)

    # MERGE DATA - allele data, alignments, resmarker info
    # Merge allele data with resmarkers, only keeping the loci we need and markers with calls
    allele_data_per_resmarker = allele_data.merge(
        processed_res_markers_info, on='Locus', how='inner')
    # Merge with alignment data
    allele_data_per_resmarker = allele_data_per_resmarker.merge(
        aligned_asv_data, left_on=['SampleID', 'Locus', 'ASV'],
        right_on=['sampleID', 'refid', 'asv'], how='left')

    # GET TABLE TO EXTRACT INFO FROM - compress to unique that need to be processed
    # Compile info for unique asv and markers
    required_columns = ['Locus', 'ASV', 'PseudoCIGAR', 'strand', 'GeneID', 'Gene', 'CodonID', 'CodonStart',
                        'hapseq', 'refseq', 'relativeCodonStart', 'CodonEnd', 'RefCodon', 'RefAA',]
    unique_asvs_per_resmarker_to_process = allele_data_per_resmarker[required_columns].drop_duplicates(
    )
    # Merge with masked coordinates
    unique_asvs_per_resmarker_to_process = unique_asvs_per_resmarker_to_process.merge(
        masked_coords_table, on='Locus', how='left')

    # THIS WAS IN A FUNCTION
    # GET UNIQUE RESMARKERS
    # Get Codon
    # This runs on unique asvs
    resmarker_per_unique_asv = unique_asvs_per_resmarker_to_process.apply(
        CodonProcessor.get_codon_info_for_marker, axis=1)
    # Add on additional columns
    resmarker_per_unique_asv['AA'] = resmarker_per_unique_asv['Codon'].apply(
        lambda x: translate(x))
    resmarker_per_unique_asv['CodonRefAlt'] = np.where(
        resmarker_per_unique_asv['RefCodon'] == resmarker_per_unique_asv['Codon'], 'REF', 'ALT')
    resmarker_per_unique_asv['AARefAlt'] = np.where(
        resmarker_per_unique_asv['RefAA'] == resmarker_per_unique_asv['AA'], 'REF', 'ALT')

    # Merge back with sample data
    resmarker_data = resmarker_per_unique_asv.merge(
        allele_data_per_resmarker[['SampleID', 'Locus', 'ASV', 'PseudoCIGAR', 'Reads', 'relativeCodonStart']], on=['Locus', 'ASV', 'PseudoCIGAR', 'relativeCodonStart'])

    # Create microhaplotype version
    mhap_table = generate_mhap_table(resmarker_data)

    # Sum reads for duplicate codons
    resmarker_data = resmarker_data.groupby(['SampleID', 'Locus', 'RefCodon', 'Codon', 'GeneID', 'Gene', 'CodonID',
                                             'relativeCodonStart', 'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt', 'strand']).Reads.sum().reset_index()

    # Ensure reads is integer
    resmarker_data['Reads'] = resmarker_data['Reads'].astype(int)

    # Sum reads for duplicate codons
    mhap_table = mhap_table.groupby(['SampleID', 'Locus', 'GeneID', 'Gene', 'MicrohapIndex',
                                    'RefMicrohap', 'Microhaplotype', 'MicrohapRefAlt']).Reads.sum().reset_index()

    # Ensure reads is integer
    mhap_table['Reads'] = mhap_table['Reads'].astype(int)

    # Sort by CodonID, Gene, and SampleID (and Locus if applicable)
    logging.debug(f"Sorting by columns")

    resmarker_data.sort_values(
        by=['SampleID', 'Locus', 'CodonID'], inplace=True)
    # Order columns
    resmarker_data = resmarker_data[['SampleID', 'Locus', 'GeneID', 'Gene', 'CodonID',
                                     'RefCodon', 'Codon', 'relativeCodonStart', 'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt', 'Reads']]

    # Collapse markers on covered by tiled amplicons
    resmarker_data_collapsed = resmarker_data.groupby(['SampleID', 'RefCodon', 'Codon', 'GeneID', 'Gene', 'CodonID',
                                                       'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt',]).Reads.sum().reset_index()
    resmarker_data_collapsed['MultipleLoci'] = resmarker_data.duplicated(
        ['SampleID', 'RefCodon', 'Codon', 'GeneID', 'Gene', 'CodonID', 'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt'], keep=False)

    # Get all mutations
    all_mutations = generate_all_mutations_table(unique_asvs_per_resmarker_to_process[[
        'Locus', 'ASV', 'PseudoCIGAR', 'strand', 'GeneID', 'Gene', 'hapseq', 'refseq']].drop_duplicates(), allele_data, ref_sequences)
    # Ensure reads is integer
    all_mutations['Reads'] = all_mutations['Reads'].astype(int)

    resmarker_data_collapsed.to_csv('resmarker_table.txt',
                                    sep='\t', index=False)
    resmarker_data.rename(
        columns={'relativeCodonStart': 'CodonStart'}, inplace=True
    )
    resmarker_data.to_csv('resmarker_table_by_locus.txt',
                          sep='\t', index=False)
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
# TODO: format grouping of the outputs
# TODO: turn these into classes
# TODO: decide on position to show : either strand and position realtive to start or just the one in the other table (or both)
# TODO: Check if masking works 0 based?
# TODO: new mutations only include unmasked ones - so a masked resmarker wouldn't be in this file
print(datetime.now() - startTime)
