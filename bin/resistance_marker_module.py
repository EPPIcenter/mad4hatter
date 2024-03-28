#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from Bio import SeqIO
from Bio.Seq import translate
import argparse
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import numpy as np
import json
import logging
from logging.handlers import RotatingFileHandler
import traceback


def parse_pseudo_cigar(string, orientation, refseq_len):
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

    logging.debug(f"Entering `parse_pseudo_cigar` with string={string}, orientation={orientation}, refseq_len={refseq_len}")

    if not isinstance(string, str):
        logging.error(f"Expected a string for `parse_pseudo_cigar` but got {type(string)}. Value: {string}")
        raise TypeError(f"Expected a string for `parse_pseudo_cigar` but got {type(string)}. Value: {string}")

    string = string.replace("=", "")  # Remove the "=" separator
    logging.debug(f"Removed '=' separator. Modified string={string}")

    transtab = str.maketrans("TACG", "ATGC")
    tuples = []

    pattern = re.compile(r'(\d+)([ID\+]?)\+?(\d+|[ACGTN]+)?')
    matches = pattern.findall(string)
    logging.debug(f"Pattern found {len(matches)} matches in the string")

    for match in matches:
        position, operation, value = match 
        position = refseq_len - int(position) if orientation == "-" else int(position) - 1 # use base 0 indexing
        logging.debug(f"Processing match: original position={match[0]}, adjusted position={position}, operation={operation}, value={value}")

        # Handle mask
        if operation == "+":  
            logging.debug("Operation is a mask. No changes applied.")
            pass # do nothing for now
        else:     
            value = value.translate(transtab) if orientation == "-" else value
            tuples.append((position, None if operation == '' else operation, value))
            logging.debug(f"Added tuple to result: {(position, None if operation == '' else operation, value)}")

    logging.debug(f"Exiting `parse_pseudo_cigar` with {len(tuples)} tuples generated")
    return tuples



def calculate_aa_changes(row, ref_sequences) -> dict:
    """
    Uses the pseudo cigar string to determine if the codon and/or amino acid matches the reference
    :param row: The currently processed row of the merged allele_data + resistance marker table
    :param ref_sequences: A dictionary of reference sequences
    :return: row 
	
    """
    logging.debug(f"------ Start of `calculate_aa_changes` for row {row['Locus']} ------")
    
    PseudoCIGAR = row['PseudoCIGAR']
    orientation = row['strand']
    logging.debug(f"Processing PseudoCIGAR: {PseudoCIGAR}, orientation: {orientation}")

    if not isinstance(PseudoCIGAR, str):
        raise TypeError(f"Expected a string for `calculate_aa_changes` but got {type(PseudoCIGAR)}. Value: {PseudoCIGAR}")

    if not isinstance(orientation, str):
        raise TypeError(f"Expected a string for `calculate_aa_changes` but got {type(orientation)}. Value: {orientation}")

    new_mutations = {}

    if PseudoCIGAR == ".":
        row['Codon'] = row['RefCodon']
        row['AA'] = translate(row['Codon'])
        row['CodonRefAlt'] = 'REF'
        row['AARefAlt'] = 'REF'
    else:
        refseq_len = len(ref_sequences[row['Locus']].seq)
        changes = parse_pseudo_cigar(PseudoCIGAR, orientation, refseq_len)
        logging.debug(f"Parsed changes: {changes}")

        codon = list(row['RefCodon'])
        # build the ASV codon using the reference codon and the changes listed in the cigar string

        # we need to eventually use the operation informatin to make informed decisions on the codon in terms of frame shifts
        for pos, op, alt in changes:
            logging.debug(f"Processing change: pos={pos}, op={op}, alt={alt}")

            # ignore masks
            if op == "+":
                continue

            # make sure that the position is within the codon
            if (pos >= row['CodonStart']) and (pos < row['CodonEnd']):
                previous_codon = codon
                index = pos - row['CodonStart']
                codon[index] = alt
                logging.debug(f"Changing codon position {index} to {alt}. Previous codon: {''.join(previous_codon)}, current codon: {''.join(codon)}")
            else:
                new_mutations[pos] = (alt, ref_sequences[row['Locus']].seq[int(pos)-1])
                logging.debug(f"Position {pos} outside of codon. Adding to new_mutations.")

        row['Codon'] = "".join(codon)
        row['AA'] = translate(row['Codon'])
        row['CodonRefAlt'] = 'ALT' if row['Codon'] != row['RefCodon'] else 'REF'
        row['AARefAlt'] = 'ALT' if row['AA'] != translate(row['RefCodon']) else 'REF'
        logging.debug(f"Final codon: {row['Codon']}, Final AA: {row['AA']}")


    # Add new mutations to the row
    row['new_mutations'] = json.dumps(new_mutations)
    logging.debug(f"New mutations added to row: {row['new_mutations']}")
    
    logging.debug(f"------ End of `calculate_aa_changes` for row {row['Locus']} ------")
    return row


def extract_info_from_gene_id(v5_string) -> tuple:
    """
    Extracts information from the gene_id column of the resistance marker table
    :param v5_string: the gene_id string found in the resistance marker table
    :return: (gene_id, gene, codon_id) 
    Example: extract_info_from_V5("PF3D7_0709000-crt-73") -> ("0709000", "crt", "73")
	
    """
    logging.debug(f"Entering `extract_info_from_V5` with v5_string={v5_string}")
    split_string = v5_string.split('-')
    gene_id = split_string[0].split('_')[-1]
    gene = split_string[1]
    codon_id = split_string[-1]
    return gene_id, gene, codon_id


def process_row(row, ref_sequences):
    logging.debug(f"Entering `process_row` with SampleID={row['SampleID']}, pseudocigar={row['PseudoCIGAR']}, Reads={row['Reads']}")
    gene_id, gene, codon_id = extract_info_from_gene_id(row['gene_id'])
    row['GeneID'] = gene_id
    row['Gene'] = gene
    row['CodonID'] = codon_id

    # Get codon and translate
    if row['strand'] == '-':
        refseq_len = len(ref_sequences[row['Locus']].seq)
        refseq_rc = ref_sequences[row['Locus']].seq.reverse_complement()
        CodonStart = refseq_len - (row['CodonStart']) - 2 # 0-based indexing
        codon_end = refseq_len - (row['CodonStart']) + 1

        row['CodonStart'] = CodonStart
        row['CodonEnd'] = codon_end

        ref_codon = str(refseq_rc[row['CodonStart'] : row['CodonEnd']])
        row['RefCodon'] = ref_codon
        row['RefAA'] = translate(ref_codon)
    else:
        CodonStart = row['CodonStart']-1
        codon_end = row['CodonStart']+2

        row['CodonStart'] = CodonStart
        row['CodonEnd'] = codon_end
        ref_codon = str(ref_sequences[row['Locus']].seq[row['CodonStart'] : row['CodonEnd']])
        row['RefCodon'] = ref_codon
        row['RefAA'] = translate(ref_codon)

    return calculate_aa_changes(row, ref_sequences)


def main(args):

    logging.debug(f"------ Start of `main` ------")
    
    # Reading allele data
    logging.debug(f"Reading allele data from: {args.allele_data_path}")
    allele_data = pd.read_csv(args.allele_data_path, sep='\t')
    logging.debug(f"Read {allele_data.shape[0]} rows from allele data")

    # Reading resistance markers
    logging.debug(f"Reading res markers info from: {args.res_markers_info_path}")
    res_markers_info = pd.read_csv(args.res_markers_info_path, sep='\t')
    logging.debug(f"Read {res_markers_info.shape[0]} rows from res markers info")

    # Renaming columns
    logging.debug("Renaming columns in res markers info")
    res_markers_info = res_markers_info.rename(columns={'codon_start': 'CodonStart', 'amplicon': 'Locus'})

    # Keep the data that we are interested in
    res_markers_info = res_markers_info[(res_markers_info['CodonStart'] > 0) & (res_markers_info['CodonStart'] < res_markers_info['ampInsert_length'])]

    # Filter allele data to only include drug resistance amplicons
    logging.debug(f"Filtering allele data based on locus suffix")
    allele_data = allele_data[allele_data['Locus'].str.endswith(('-1B', '-2'))]
    logging.debug(f"Filtered allele data to {allele_data.shape[0]} rows")

    # Join the allele table and resistance marker table on the locus
    logging.debug("Joining allele table and resistance marker table on the locus")
    allele_data = res_markers_info.set_index('Locus').join(allele_data.set_index('Locus'), how='left').reset_index()
    logging.debug(f"Joined data has {allele_data.shape[0]} rows")

    # Filter any rows that have `NaN` in the SampleID column 
    allele_data = allele_data.dropna(subset=['PseudoCIGAR'])

    # Ensure reads is integer
    allele_data['Reads'] = allele_data['Reads'].astype(int)

    # Read in the reference sequences - this will be used for the reference codons
    logging.debug("Reading in the reference sequences")
    ref_sequences = SeqIO.to_dict(SeqIO.parse(args.refseq_path, 'fasta'))
    logging.debug(f"Loaded {len(ref_sequences)} reference sequences")

    # Create function to run on all rows of the joined allele data + resistance marker table.
    logging.debug("Setting up the parallel function")
    process_row_partial = partial(process_row, ref_sequences=ref_sequences)

    # Run in parallel
    logging.debug(f"Processing rows in parallel using {args.n_cores} cores")

    results = []

    with ProcessPoolExecutor(max_workers=args.n_cores) as executor:
        # Create a dictionary that maps futures to their corresponding rows for debugging and error handling
        future_to_row = {executor.submit(process_row_partial, row[1]): row for row in allele_data.iterrows()}

        for future in as_completed(future_to_row):
            row_data = future_to_row[future]
            try:
                result = future.result()
                results.append(result)
                locus_val = row_data[1]['Locus']
                sample_id = row_data[1]['SampleID']
                reads = row_data[1]['Reads']
                gene_id = row_data[1]['gene_id']
                logging.debug(f"Successfully processed row with SampleID: {sample_id}, Locus: {locus_val}, Reads: {reads}, GeneID: {gene_id}")
            except Exception as e:
                logging.error(f"Error processing row with locus: {row_data[1]['Locus']}. Error: {e}")

    
    logging.debug(f"Total rows processed: {len(results)}")
    # Create a table that identifies whether there are dna and/or codon difference in the ASVs
    # at the positions specified in the resistance marker table
    df_results = pd.DataFrame(results)
    df_results = df_results.rename(columns={'reads': 'Reads', 'sampleID': 'SampleID', "pseudo_cigar": "PseudoCIGAR"})
    df_results['CodonID'] = df_results['CodonID'].astype(int)

    sort_columns, sort_order = ['SampleID', 'chr', 'CodonID'], [True, True, True]
    select_columns = ['SampleID', 'chr', 'Locus', 'GeneID', 'Gene', 'CodonID', 'RefCodon', 'Codon', 'CodonStart', 'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt', 'Reads', 'PseudoCIGAR', 'new_mutations']

    # Select the columns we want
    logging.debug(f"Selecting columns: {select_columns}")
    df_resmarker = df_results[select_columns]
    logging.debug(f"Selected {df_results.shape[0]} rows")
    logging.debug(f"Columns: {df_results.columns}")

    # Summarize reads
    group_by_columns = ['SampleID', 'Locus', 'chr', 'GeneID', 'Gene', 'CodonID', 'RefCodon', 'Codon', 'CodonStart', 'CodonRefAlt', 'RefAA', 'AA', 'AARefAlt']

    logging.debug(f"Group by columns for resmarker table: {group_by_columns}")

    df_resmarker = df_resmarker.groupby(group_by_columns).agg({
        'Reads': 'sum'
    }).reset_index()

    # Sort by CodonID, Gene, and SampleID (and Locus if applicable)
    logging.debug(f"Sorting by columns: {sort_columns}, order: {sort_order}")
    df_resmarker.sort_values(by=sort_columns, ascending=sort_order, inplace=True)

    # Drop unnecessary columns
    drop_columns = ['chr']
    logging.debug(f"Dropping columns: {drop_columns}")
    df_resmarker.drop(drop_columns, axis=1, inplace=True)

    # Output resmarker table
    df_resmarker.to_csv('resmarker_table.txt', sep='\t', index=False)

    def create_microhap_lambda(x):
        # Obtain sorted codons and their corresponding amino acids
        sorted_codon_ids = x['CodonID'].sort_values()
        amino_acids_by_codon = x.set_index('CodonID').loc[sorted_codon_ids]['AA']

        # Check for length mismatch
        if len(sorted_codon_ids) != len(amino_acids_by_codon):
            logging.error(f"Length mismatch detected! Codon IDs: {sorted_codon_ids.tolist()}, Amino Acids: {amino_acids_by_codon.tolist()}")

        return pd.Series({
            'MicrohapIndex': '/'.join(map(str, sorted_codon_ids)),
            'Microhaplotype': '/'.join(amino_acids_by_codon),
            'RefMicrohap': '/'.join(x.set_index('CodonID').loc[sorted_codon_ids]['RefAA']),
        })

    group_by_columns = ['SampleID', 'Locus', 'chr', 'GeneID', 'Gene', 'PseudoCIGAR', 'Reads']

    df_microhap = df_results.groupby(group_by_columns).apply(create_microhap_lambda).reset_index()

    # Create MicrohapRefAlt column
    df_microhap['MicrohapRefAlt'] = np.where(df_microhap['Microhaplotype'] == df_microhap['RefMicrohap'], 'REF', 'ALT')

    # Select columns and rename them
    microhap_select_columns = ['SampleID', 'Locus', 'chr', 'GeneID', 'Gene', 'MicrohapIndex', 'RefMicrohap', 'Microhaplotype', 'MicrohapRefAlt', 'Reads']
    df_microhap = df_microhap[microhap_select_columns]

    # Summarize reads
    microhap_groupby_columns = ['SampleID', 'Locus', 'chr', 'GeneID', 'Gene', 'MicrohapIndex', 'RefMicrohap', 'Microhaplotype', 'MicrohapRefAlt']

    df_microhap_collapsed = df_microhap.groupby(microhap_groupby_columns).agg({
        'Reads': 'sum'
    }).reset_index()

    # Function to extract the first integer from MicrohapIndex for sorting
    def extract_first_position(microhap_index):
        # Split the string, convert to integers, and return the first element
        positions = list(map(int, microhap_index.split('/')))
        return min(positions)

    df_microhap_collapsed['SortPosition'] = df_microhap_collapsed['MicrohapIndex'].apply(extract_first_position)

    # Sort by MicrohapIndex, Gene, and SampleID
    df_microhap_collapsed_sort_columns, df_microhap_collapsed_sort_order = ['SampleID', 'chr', 'SortPosition'], [True, True, True]
    df_microhap_collapsed = df_microhap_collapsed.sort_values(by=df_microhap_collapsed_sort_columns, ascending=df_microhap_collapsed_sort_order)

    # Remove the sort position column
    df_microhap_collapsed.drop(['SortPosition', 'chr'], axis=1, inplace=True)

    # Output microhaplotype table
    df_microhap_collapsed.to_csv('resmarker_microhap_table.txt', sep='\t', index=False)

    # Prepare the data for the new mutations output
    mutation_list = []
    for _, row in df_results.iterrows():
        new_mutations = json.loads(row['new_mutations'])
        if len(new_mutations) > 0:
            for pos, (alt, ref) in new_mutations.items():
                new_row = {
                    'SampleID': row['SampleID'],
                    'Locus': row['Locus'],
                    'GeneID': row['GeneID'],
                    'Gene': row['Gene'],
                    'CodonID': row['CodonID'],
                    'Position': pos,
                    'Alt': alt,
                    'Ref': ref,
                    'Reads': row['Reads']
                }

                mutation_list.append(new_row)

    df_new_mutations = pd.DataFrame(mutation_list)
    df_new_mutations.to_csv('resmarker_new_mutations.txt', sep='\t', index=False)


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="Process some input files.")

        parser.add_argument("--allele_data_path", required=True, help="Path to allele_data.txt")
        parser.add_argument("--res_markers_info_path", required=True, help="Path to resistance marker table")
        parser.add_argument("--refseq_path", required=True, help="Path to reference sequences [fasta]")
        parser.add_argument("--n-cores", type=int, default=1, help="Number of cores to use")
        parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Set the logging level")
        parser.add_argument("--log-file", type=str, help="Write logs to a file instead of the console")
        parser.add_argument("--log-max-size", type=int, default=5, help="Maximum log file size in MB")
        parser.add_argument("--log-backups", type=int, default=1, help="Number of log files to keep")

        args = parser.parse_args()

        numeric_level = getattr(logging, args.log_level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f"Invalid log level: {args.log_level}")

        if args.log_file:
            # Use rotating log files
            file_handler = RotatingFileHandler(args.log_file, maxBytes=args.log_max_size * 1024 * 1024, backupCount=args.log_backups)
            logging.basicConfig(level=numeric_level, format="%(asctime)s - %(levelname)s - %(message)s", handlers=[file_handler])
        else:
            logging.basicConfig(level=numeric_level, format="%(asctime)s - %(levelname)s - %(message)s")

        logging.info("Program start.")
        main(args)
    except Exception as e:
        logging.error(f"An error of type {type(e).__name__} occurred. Message: {e}")
        logging.error(traceback.format_exc())  # Log the full traceback
        raise e
    finally:
        logging.info("Program complete.")
        logging.shutdown()

