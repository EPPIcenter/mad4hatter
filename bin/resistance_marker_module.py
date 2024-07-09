#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from Bio import SeqIO
from Bio.Seq import translate, Seq
import argparse
import re
import numpy as np
import logging
from logging.handlers import RotatingFileHandler
import traceback


def read_allele_data(allele_data_path, aligned_asv_table_path):
    """Reads and merges allele data with aligned ASV data."""
    logging.info(
        f"Reading allele data from {allele_data_path} and aligned ASV data from {aligned_asv_table_path}.")
    try:
        allele_data = pd.read_csv(allele_data_path, sep='\t')
        aligned_asv_data = pd.read_csv(aligned_asv_table_path, sep='\t')
        merged_data = allele_data.merge(
            aligned_asv_data, left_on=['SampleID', 'Locus', 'ASV'],
            right_on=['sampleID', 'refid', 'asv'], how='left'
        )
        logging.info("Allele data and aligned ASV data merged successfully.")
        return merged_data
    except Exception as e:
        logging.error(f"Failed to read or merge data: {e}")
        raise


class ProcessMarkerInfo:
    def __init__(self, ref_sequences):
        self.ref_sequences = ref_sequences

    def extract_marker_info(self, row):
        """Adds marker start position relative to strand and reference to a row of marker information."""
        logging.debug(
            f"Extracting marker information for row with Locus {row['Locus']}")

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
            'CodonStart': codon_start,
            'relativeCodonStart': relative_codon_start,
            'CodonEnd': codon_end,
            'RefCodon': ref_codon
        }

    def extract_position_and_reference_for_markers(self, res_marker_table):
        """Extracts position and reference information for markers."""
        logging.info(
            f"Extracting position and reference information for {len(res_marker_table)} markers.")
        marker_info = res_marker_table.apply(self.extract_marker_info, axis=1)
        res_marker_table = pd.DataFrame(marker_info.tolist())
        res_marker_table['RefAA'] = res_marker_table['RefCodon'].apply(
            lambda x: 'X' if not translate(x) else translate(x)
        )
        logging.info(
            "Position and reference information extracted successfully.")
        return res_marker_table

    def check_if_marker_masked(self, resmarker_table, mask_coords_dict):
        # Generate 1-based codon coordinates relative to forward strand
        resmarker_table['CodonCoordinates'] = resmarker_table.apply(
            lambda row: set(range(row.CodonStart - 1, row.CodonStart + 2)), axis=1)

        # Function to check if codon is masked
        def is_codon_masked(row):
            mask_coords = mask_coords_dict.get(row.Locus, set())
            return not mask_coords.isdisjoint(row.CodonCoordinates)

        resmarker_table['CodonMasked'] = resmarker_table.apply(
            is_codon_masked, axis=1)
        resmarker_table.drop(columns=['CodonCoordinates'], inplace=True)
        return resmarker_table


class MaskCoordinatesExtractor:
    @staticmethod
    def mask_str_to_coordinates(mask_string):
        """Converts a mask string to a list of 1-based masked coordinates."""
        start, length = map(int, mask_string.strip('N').split('+'))
        start = start
        end = start + length
        return list(range(start, end))

    @staticmethod
    def pseudocigar_to_mask_coordinates(pseudo_cigar):
        """Converts a pseudo CIGAR string to 1-based mask coordinates."""
        logging.debug(
            f"Extracting masking coordinates from cigar {pseudo_cigar}.")
        masking_pattern = r'\d+\+\d+N'
        masks = re.findall(masking_pattern, pseudo_cigar)
        mask_coordinates = []
        for mask in masks:
            mask_coordinates.extend(
                MaskCoordinatesExtractor.mask_str_to_coordinates(mask))
        return set(mask_coordinates)

    @staticmethod
    def create_masked_coordinate_dict(allele_data):
        """Creates a table of 1-based masked coordinates from allele data."""
        logging.info(
            f"Creating masked coordinate table for {allele_data.Locus.nunique()} loci.")
        cigar_per_locus = allele_data.groupby('Locus', as_index=False).first()
        cigar_per_locus['masked_coords'] = cigar_per_locus['PseudoCIGAR'].apply(
            MaskCoordinatesExtractor.pseudocigar_to_mask_coordinates)
        # Create the dictionary
        masked_coords_dict = cigar_per_locus.set_index(
            'Locus')['masked_coords'].to_dict()
        logging.info("Masked coordinate dictionary created successfully.")
        print(masked_coords_dict)
        return masked_coords_dict


class CodonProcessor:
    @staticmethod
    def set_codon_with_indel_as_undetermined(codon, ref_codon):
        """Sets codons with indels as undetermined."""
        return 'X' if '-' in codon or '-' in ref_codon else codon

    @staticmethod
    def extract_codon_from_asv(reference, asv, mask_coordinates, start_position):
        """Extracts codon from ASV sequence with respect to strand and indels."""
        follows_indel = False
        if '-' not in reference and '-' not in asv:
            codon = asv[start_position:start_position+3]
        else:
            index = 0
            while index < start_position:
                ref_base = reference[index]
                seq_base = asv[index]
                if ref_base == '-':
                    start_position += 1
                    if index + 1 not in mask_coordinates:
                        follows_indel = True
                if seq_base == '-':
                    if index + 1 not in mask_coordinates:
                        follows_indel = True
                index += 1
            codon = asv[start_position:start_position+3]
            aligned_ref_codon = reference[start_position:start_position+3]
            codon = CodonProcessor.set_codon_with_indel_as_undetermined(
                codon, aligned_ref_codon)
        return codon, follows_indel

    @staticmethod
    def get_codon_info_for_marker(marker_info, masking_info):
        """Gets codon information for a marker."""
        logging.debug(
            f"Getting codon information for marker at position {marker_info['relativeCodonStart']} to {marker_info['CodonEnd']} from {marker_info['hapseq']}")
        reference = marker_info['refseq']
        sequence = marker_info['hapseq']
        start_position = marker_info['relativeCodonStart']
        pseudo_cigar = marker_info['PseudoCIGAR']
        ref_codon = marker_info['RefCodon']
        mask_coordinates = masking_info.get(marker_info.Locus, set())
        strand = marker_info['strand']
        end_position = marker_info['CodonEnd']

        if pseudo_cigar == '.':
            codon = ref_codon
            follows_indel = False
        else:
            if strand == '-':
                sequence = Seq(sequence).reverse_complement()
                reference = Seq(reference).reverse_complement()
            # Get codon from asv
            codon, follows_indel = CodonProcessor.extract_codon_from_asv(
                reference, sequence, mask_coordinates, start_position
            )
        marker_with_codon_info = (pd.concat([marker_info, pd.Series(
            {'Codon': codon, 'FollowsIndel': follows_indel})]))
        return marker_with_codon_info


def pseudo_cigar_to_mutations(pseudo_cigar):
    """Parse a pseudo CIGAR string into a list of tuples (position, operation, value)."""
    logging.debug(f"Entering `parse_pseudo_cigar` with string={pseudo_cigar}")

    if not isinstance(pseudo_cigar, str):
        logging.error(
            f"Expected a string for `parse_pseudo_cigar` but got {type(pseudo_cigar)}. Value: {pseudo_cigar}")
        raise TypeError(
            f"Expected a string for `parse_pseudo_cigar` but got {type(pseudo_cigar)}. Value: {pseudo_cigar}")

    pseudo_cigar = pseudo_cigar.replace("=", "")  # Remove the "=" separator
    logging.debug(f"Removed '=' separator. Modified string={pseudo_cigar}")

    tuples = []

    pattern = re.compile(r'(\d+)([ID\+]?)\+?(\d+|[ACGTN]+)?')
    matches = pattern.findall(pseudo_cigar)
    logging.debug(f"Pattern found {len(matches)} matches in the string")

    for match in matches:
        position, operation, value = match
        logging.debug(
            f"Processing match: original position={match[0]}, adjusted position={position}, operation={operation}, value={value}")

        # Handle mask
        if operation == "+":
            logging.debug("Operation is a mask. No changes applied.")
            pass  # do nothing for now
        else:
            tuples.append((position, None if operation ==
                          '' else operation, value))
            logging.debug(
                f"Added tuple to result: {(position, None if operation == '' else operation, value)}")

    logging.debug(
        f"Exiting `parse_pseudo_cigar` with {len(tuples)} tuples generated")
    return tuples


def extract_mutations_from_unique_pseudo_cigar(df, ref_sequences):
    """Extracts mutations from pseudo CIGAR strings in dataframe."""
    logging.info("Extracting mutations from unique pseudo CIGAR strings.")
    results = []
    transtab = str.maketrans("TACG", "ATGC")
    for index, row in df.iterrows():
        pseudo_cigar = row['PseudoCIGAR']
        locus = row['Locus']
        strand = row['strand']

        if pseudo_cigar == ".":
            continue

        # returns 1-based position
        changes = pseudo_cigar_to_mutations(pseudo_cigar)

        ref_seq = ref_sequences[locus].seq
        for pos, op, alt in changes:
            pos = int(pos)
            ref = ref_seq[pos-1]
            if strand == '-':
                alt = alt.translate(transtab)
                ref = ref.translate(transtab)
            if op == 'D':
                ref = alt
                alt = '-'
            elif op == 'I':
                ref = '-'
            results.append({
                'GeneID': row.GeneID,
                'Gene': row.Gene,
                'Locus': locus,
                'PseudoCIGAR': pseudo_cigar,
                'LocusPosition': pos,
                'Alt': alt,
                'Ref': ref,
            })
    logging.info("Mutations extracted successfully.")
    return pd.DataFrame(results)


def assign_ref_alt(df, ref_colunn, data_column):
    return np.where(
        df[data_column] == 'X',
        'X',
        np.where(
            df[ref_colunn] == df[data_column], 'REF', 'ALT'
        )
    )


class ResmarkerTableGenerator:
    @staticmethod
    def generate_resmarker_table_for_unique_asv(unique_asvs_per_resmarker_to_process, masked_coords_table):
        """Build table of the resmarkers in each unique asv."""
        logging.info("Extracting resmarker info from unique ASVs.")
        resmarker_per_unique_asv = unique_asvs_per_resmarker_to_process.apply(
            CodonProcessor.get_codon_info_for_marker, masking_info=masked_coords_table,  axis=1)
        # Add on additional columns
        resmarker_per_unique_asv['AA'] = resmarker_per_unique_asv['Codon'].apply(
            lambda x: 'X' if not translate(x) else translate(x)
        )
        resmarker_per_unique_asv['CodonRefAlt'] = assign_ref_alt(
            resmarker_per_unique_asv, 'RefCodon', 'Codon')

        resmarker_per_unique_asv['AARefAlt'] = assign_ref_alt(
            resmarker_per_unique_asv, 'RefAA', 'AA')
        return resmarker_per_unique_asv

    @staticmethod
    def assign_resmarkers_to_samples(sample_allele_data, resmarker_for_asv):
        """Merge resmarker information from asv to samples with that asv."""
        logging.info("Merging resmarker info with samples.")
        # Merge back with sample data
        resmarker_data = resmarker_for_asv.merge(
            sample_allele_data[['SampleID', 'Locus', 'ASV', 'PseudoCIGAR', 'Reads', 'relativeCodonStart']], on=['Locus', 'ASV', 'PseudoCIGAR', 'relativeCodonStart'])
        # Sum reads for duplicate codons
        resmarker_data = resmarker_data.groupby(['SampleID', 'GeneID', 'Gene', 'Locus', 'CodonID', 'RefCodon',
                                                'Codon', 'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt', 'FollowsIndel', 'CodonMasked']).Reads.sum().reset_index()
        resmarker_data['Reads'] = resmarker_data['Reads'].astype(int)
        # Sort by CodonID, Gene, and SampleID (and Locus if applicable)
        logging.debug(f"Sorting by columns")
        resmarker_data.sort_values(
            by=['SampleID', 'Locus', 'CodonID'], inplace=True)
        return resmarker_data

    @staticmethod
    def collapse_resmarkers_on_tiled_amplicons(resmarker_table):
        """For markers present on multiple loci sum reads and flag from multiple loci"""
        logging.info("Collapsing tiled markers.")
        resmarker_table_by_locus = resmarker_table.copy()
        columns_to_collapse = ['SampleID', 'GeneID', 'Gene', 'CodonID', 'RefCodon',
                               'Codon', 'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt', 'FollowsIndel']
        resmarker_table_by_locus['MultipleLoci'] = resmarker_table_by_locus.duplicated(
            columns_to_collapse, keep=False)
        resmarker_data_collapsed = resmarker_table_by_locus.groupby(columns_to_collapse+['MultipleLoci']).agg({
            'CodonMasked': 'any',
            'Reads': 'sum'
        }).reset_index()
        # reorder columns
        resmarker_data_collapsed = resmarker_data_collapsed[['SampleID', 'GeneID', 'Gene', 'CodonID', 'RefCodon', 'Codon',
                                                             'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt', 'FollowsIndel', 'CodonMasked', 'MultipleLoci', 'Reads']]
        return resmarker_data_collapsed

    @staticmethod
    def generate_mhap_table_for_unique_asv(resmarker_table):
        """Generate a table of microhaplotypes for markers on locus for each unique asv."""
        logging.info("Compiling microhaplotype of resmarkers for unique ASVs.")
        resmarker_table = resmarker_table.sort_values(by='relativeCodonStart')
        mhap_table = resmarker_table.groupby(['ASV', 'Locus', 'GeneID', 'Gene', 'PseudoCIGAR']).agg(
            MicrohaplotypeCodonIDs=('CodonID',
                                    lambda x: '/'.join(map(str, x))),
            RefMicrohap=('RefAA', lambda x: '/'.join(map(str, x))),
            Microhaplotype=('AA', lambda x: '/'.join(map(str, x)))
        ).reset_index()
        mhap_table['MicrohapRefAlt'] = assign_ref_alt(
            mhap_table, 'RefMicrohap', 'Microhaplotype')
        return mhap_table

    @staticmethod
    def assign_resmarker_mhaps_to_samples(sample_allele_data, mhap_for_asv):
        """Merging microhaplotype information for markers back to samples."""
        logging.info("Merging microhaplotype information with samples.")
        mhap_table = mhap_for_asv.merge(
            sample_allele_data[['SampleID', 'Locus', 'ASV', 'PseudoCIGAR', 'Reads']], on=['Locus', 'ASV', 'PseudoCIGAR'])
        mhap_table = mhap_table.groupby(['SampleID',  'GeneID',  'Gene', 'Locus', 'MicrohaplotypeCodonIDs',
                                        'RefMicrohap', 'Microhaplotype', 'MicrohapRefAlt']).Reads.sum().reset_index()
        mhap_table['Reads'] = mhap_table['Reads'].astype(int)
        return mhap_table

    @staticmethod
    def generate_all_mutations_table(unique_pseudo_cigars, allele_data, ref_sequences):
        """Generate table of all differences noted in the pseudo cigar."""
        logging.info(
            "Generating table of all mutation from unique pseudo cigars.")
        # Get all mutations
        mutations_df = extract_mutations_from_unique_pseudo_cigar(
            unique_pseudo_cigars, ref_sequences)
        all_mutations = allele_data.merge(mutations_df, on=['Locus', 'PseudoCIGAR'])[
            ['SampleID', 'Locus', 'GeneID', 'Gene', 'PseudoCIGAR', 'LocusPosition', 'Alt', 'Ref', 'Reads']]
        all_mutations = all_mutations.groupby(
            ['SampleID', 'GeneID', 'Gene', 'Locus', 'LocusPosition', 'Alt', 'Ref']).Reads.sum().reset_index()
        # Ensure reads is integer
        all_mutations['Reads'] = all_mutations['Reads'].astype(int)
        all_mutations.sort_values(
            by=['SampleID', 'Locus', 'LocusPosition'], inplace=True)
        return all_mutations

    @staticmethod
    def generate_resmarker_tables(unique_asvs_per_resmarker_df, allele_data_per_resmarker, allele_data, ref_sequences, masked_coords):
        """Generate tables of resmarkers, mhaps, and all mutations."""
        logging.info("Generating resistance marker tables.")
        # resmarker table
        resmarker_per_unique_asv = ResmarkerTableGenerator.generate_resmarker_table_for_unique_asv(
            unique_asvs_per_resmarker_df, masked_coords)
        resmarker_data = ResmarkerTableGenerator.assign_resmarkers_to_samples(
            allele_data_per_resmarker, resmarker_per_unique_asv)

        # resmarker table with tiled markers collapsed
        resmarker_data_collapsed_tiled = ResmarkerTableGenerator.collapse_resmarkers_on_tiled_amplicons(
            resmarker_data)

        # mhap table
        mhap_per_unique_asv = ResmarkerTableGenerator.generate_mhap_table_for_unique_asv(
            resmarker_per_unique_asv)
        mhap_table = ResmarkerTableGenerator.assign_resmarker_mhaps_to_samples(
            allele_data_per_resmarker, mhap_per_unique_asv
        )

        # all mutations table
        unique_pseudo_cigars = unique_asvs_per_resmarker_df[[
            'Locus', 'PseudoCIGAR', 'strand', 'GeneID', 'Gene', 'hapseq', 'refseq']].drop_duplicates()
        all_mutations = ResmarkerTableGenerator.generate_all_mutations_table(
            unique_pseudo_cigars, allele_data, ref_sequences)
        return resmarker_data, mhap_table, resmarker_data_collapsed_tiled, all_mutations


def main(args):
    """Main function to process data."""
    logging.debug(f"------ Start of `main` ------")

    # READ DATA
    allele_data = read_allele_data(
        args.allele_data_path, args.aligned_asv_table_path)
    ref_sequences = SeqIO.to_dict(SeqIO.parse(args.refseq_path, 'fasta'))
    res_markers_info = pd.read_csv(args.res_markers_info_path, dtype={
                                   'GeneID': str}, sep='\t')

    # Add on reference codon and amino acid and position relative to strand
    marker_info_extractor = ProcessMarkerInfo(ref_sequences)
    res_markers_info = marker_info_extractor.extract_position_and_reference_for_markers(
        res_markers_info)

    # Add flag if the marker is masked
    masked_coords_dict = MaskCoordinatesExtractor.create_masked_coordinate_dict(
        allele_data)
    marker_info_extractor.check_if_marker_masked(
        res_markers_info, masked_coords_dict)
    allele_data_per_resmarker = allele_data.merge(
        res_markers_info, on='Locus')

    # Get unique asv information to process
    required_columns = ['Locus', 'ASV', 'PseudoCIGAR', 'strand', 'GeneID', 'Gene', 'CodonID',
                        'hapseq', 'refseq', 'CodonStart', 'relativeCodonStart', 'CodonEnd', 'RefCodon', 'RefAA', 'CodonMasked']
    unique_asvs_per_resmarker_df = allele_data_per_resmarker[required_columns].drop_duplicates(
    )

    # Generate tables
    resmarker_data, mhap_table, resmarker_data_collapsed_tiled, all_mutations = ResmarkerTableGenerator.generate_resmarker_tables(
        unique_asvs_per_resmarker_df, allele_data_per_resmarker, allele_data, ref_sequences, masked_coords_dict)

    logging.info(f"Writing output data to files.")
    resmarker_data_collapsed_tiled.to_csv('resmarker_table.txt',
                                          sep='\t', index=False)
    resmarker_data.to_csv('resmarker_table_by_locus.txt',
                          sep='\t', index=False)
    mhap_table.to_csv('resmarker_microhaplotype_table.txt',
                      sep='\t', index=False)
    all_mutations.to_csv('all_mutations_table.txt', sep='\t', index=False)
    logging.info(f"Finished writing outputs.")


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
        logging.error(
            f"An error of type {type(e).__name__} occurred. Message: {e}")
        logging.error(traceback.format_exc())  # Log the full traceback
        raise e
    finally:
        logging.info("Program complete.")
        logging.shutdown()
