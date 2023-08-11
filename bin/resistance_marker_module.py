#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from Bio import SeqIO
from Bio.Seq import translate
import argparse
import re
from concurrent.futures import ProcessPoolExecutor
from functools import partial
import numpy as np
import json


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

    string = string.replace("=", "")  # Remove the "=" separator
    transtab = str.maketrans("TACG", "ATGC")
    tuples = []

    pattern = re.compile(r'(\d+)([ID\+]?)\+?(\d+|[ACGTN]+)?')
    matches = pattern.findall(string)
    for match in matches:
        position, operation, value = match 
        position = refseq_len - int(position) - 1 if orientation == "-" else int(position) - 1 # use base 0 indexing

        # Handle mask
        if operation == "+":  
            mask_positions = [(position + i, "N", None) for i in range(int(value))]
            tuples.extend(mask_positions)
        else:     
            value = value.translate(transtab) if orientation == "-" else value
            tuples.append((position, None if operation == '' else operation, value))

    # Added masking separately
    pattern = re.compile(r'(\d+)N\+(\d+)')
    matches = pattern.findall(string)

    return tuples


def calculate_aa_changes(row, ref_sequences) -> dict:
    """
    Uses the pseudo cigar string to determine if the codon and/or amino acid matches the reference

    :param row: The currently processed row of the merged allele_data + resistance marker table
    :param ref_sequences: A dictionary of reference sequences
    :return: row 
	
    """
    
    pseudo_cigar = row['pseudo_cigar']
    orientation = row['V4']
    new_mutations = {}

    if pseudo_cigar == ".":
        row['Codon'] = row['Reference_Codon']
        row['AA'] = translate(row['Codon'])
        row['Codon_Ref/Alt'] = 'REF'
        row['AA_Ref/Alt'] = 'REF'
    else:
        refseq_len = len(ref_sequences[row['amplicon']].seq)
        changes = parse_pseudo_cigar(pseudo_cigar, orientation, refseq_len)
        codon = list(row['Reference_Codon'])
        # build the ASV codon using the reference codon and the changes listed in the cigar string

        # we need to eventually use the operation informatin to make informed decisions on the codon in terms of frame shifts
        for pos, op, alt in changes:

            # ignore masks
            if op != "+":
                continue

            # make sure that the position is within the codon
            if (pos >= row['Codon_Start']) and (pos < row['Codon_End']):
                index = pos - row['Codon_Start']
                codon[index] = alt # replace the reference base with the alternate base
                # note: we expect substitutions in these regions
            else: 
                # if the position is outside of the codon, then we need to add the mutation to the new_mutations dictionary
                new_mutations[pos] = (alt, ref_sequences[row['amplicon']].seq[int(pos)-1])
            

        row['Codon'] = "".join(codon) # Collapse the bases into a 3 character string
        row['AA'] = translate(row['Codon']) # Use the Bio package to translate the codon to an amino acid
        row['Codon_Ref/Alt'] = 'ALT' if row['Codon'] != row['Reference_Codon'] else 'REF'
        row['AA_Ref/Alt'] = 'ALT' if row['AA'] != translate(row['Reference_Codon']) else 'REF'


    # Add new mutations to the row
    row['new_mutations'] = json.dumps(new_mutations)
    
    return row


def extract_info_from_V5(v5_string) -> tuple:
    """
    Extracts information from the V5 column of the resistance marker table

    :param v5_string: the V5 string found in the resistance marker table
    :return: (gene_id, gene, codon_id) 

    Example: extract_info_from_V5("PF3D7_0709000-crt-73") -> ("0709000", "crt", "73")
	
    """
    split_string = v5_string.split('-')
    gene_id = split_string[0].split('_')[-1]
    gene = split_string[1]
    codon_id = split_string[-1]
    return gene_id, gene, codon_id


def process_row(row, ref_sequences):
    gene_id, gene, codon_id = extract_info_from_V5(row['V5'])
    row['Gene_ID'] = gene_id
    row['Gene'] = gene
    row['Codon_ID'] = codon_id

    # Get codon and translate
    if row['V4'] == '-':
        refseq_len = len(ref_sequences[row['amplicon']].seq)
        refseq_rc = ref_sequences[row['amplicon']].seq.reverse_complement()
        codon_start = refseq_len - (row['Codon_Start']) - 2 # 0-based indexing
        codon_end = refseq_len - (row['Codon_Start']) + 1

        row['Codon_Start'] = codon_start
        row['Codon_End'] = codon_end

        ref_codon = str(refseq_rc[row['Codon_Start'] : row['Codon_End']])
        row['Reference_Codon'] = ref_codon
        row['Reference_AA'] = translate(ref_codon)
    else:
        codon_start = row['Codon_Start']-1
        codon_end = row['Codon_Start']+2

        row['Codon_Start'] = codon_start
        row['Codon_End'] = codon_end
        ref_codon = str(ref_sequences[row['amplicon']].seq[row['Codon_Start'] : row['Codon_End']])
        row['Reference_Codon'] = ref_codon
        row['Reference_AA'] = translate(ref_codon)

    return calculate_aa_changes(row, ref_sequences)


def main(args):
    allele_data = pd.read_csv(args.allele_data_path, sep='\t')
    res_markers_info = pd.read_csv(args.res_markers_info_path, sep='\t')

    # Keep the data that we are interested in
    res_markers_info = res_markers_info[(res_markers_info['Codon_Start'] > 0) & (res_markers_info['Codon_Start'] < res_markers_info['ampInsert_length'])]
    res_markers_info = res_markers_info.drop_duplicates(subset='V5', keep='first')

    # Filter allele data to only include drug resistance amplicons
    allele_data = allele_data[allele_data['locus'].str.endswith(('-1B', '-2'))]

    # Join the allele table and resistance marker table on the locus
    allele_data = res_markers_info.set_index('amplicon').join(allele_data.set_index('locus'), how='left').reset_index()

    # Read in the reference sequences - this will be used for the reference codons
    ref_sequences = SeqIO.to_dict(SeqIO.parse(args.refseq_path, 'fasta'))

    # Create function to run on all rows of the joined allele data + resistance marker table.
    process_row_partial = partial(process_row, ref_sequences=ref_sequences)

    # Run in parallel
    with ProcessPoolExecutor(max_workers=args.n_cores) as executor:
        results = list(executor.map(process_row_partial, [row for _, row in allele_data.iterrows()]))

    # Create a table that identifies whether there are dna and/or codon difference in the ASVs
    # at the positions specified in the resistance marker table
    df_results = pd.DataFrame(results)
    df_results = df_results.sort_values(['Codon_ID', 'Gene', 'sampleID'] ,ascending=[True, True, True])
    df_results = df_results[['sampleID', 'Gene_ID', 'Gene', 'Codon_ID', 'Reference_Codon', 'Codon', 'Codon_Start', 'Codon_Ref/Alt', 'Reference_AA', 'AA', 'AA_Ref/Alt', 'reads', 'amplicon', 'pseudo_cigar', 'new_mutations']]
    df_results = df_results.rename(columns={'reads': 'Reads'})
    df_results.drop(['amplicon', 'pseudo_cigar', 'new_mutations'], axis=1).to_csv('resmarker_table.txt', sep='\t', index=False)


    # Group data and create microhaplotypes
    df_microhap = df_results.groupby(['sampleID', 'Gene_ID', 'Gene', 'pseudo_cigar', 'Reads']).apply(
        lambda x: pd.Series({
            'Microhaplotype_Index': '/'.join(map(str, x['Codon_ID'].sort_values())),
            'Microhaplotype': '/'.join(x.set_index('Codon_ID').loc[x['Codon_ID'].sort_values()]['AA']),
            'Reference_Microhaplotype': '/'.join(x.set_index('Codon_ID').loc[x['Codon_ID'].sort_values()]['Reference_AA']),
        })
    ).reset_index()

    # Create Microhaplotype_Ref/Alt column
    df_microhap['Microhaplotype_Ref/Alt'] = np.where(df_microhap['Microhaplotype'] == df_microhap['Reference_Microhaplotype'], 'REF', 'ALT')

    # Select columns and rename them
    df_microhap = df_microhap[['sampleID', 'Gene_ID', 'Gene', 'Microhaplotype_Index', 'Reference_Microhaplotype', 'Microhaplotype', 'Microhaplotype_Ref/Alt', 'Reads']]
    df_microhap = df_microhap.sort_values(['Microhaplotype_Index', 'Gene', 'sampleID'], ascending=[True, True, True])
    
    # Summarize reads
    df_microhap_collapsed = df_microhap.groupby(['sampleID', 'Gene_ID', 'Gene', 'Microhaplotype_Index', 'Reference_Microhaplotype', 'Microhaplotype', 'Microhaplotype_Ref/Alt']).agg({
        'Reads': 'sum'
    }).reset_index()

    # Output microhaplotype table
    df_microhap_collapsed.to_csv('resmarker_microhap_table.txt', sep='\t', index=False)

    # Create New Mutations table (could parellize)
    mutation_list = []
    for _, row in df_results.iterrows():
        new_mutations = json.loads(row['new_mutations'])
        if len(new_mutations) > 0:
            for pos, (alt, ref) in new_mutations.items():
                new_row = {
                    'sampleID': row['sampleID'],
                    'Gene_ID': row['Gene_ID'],
                    'Gene': row['Gene'],
                    'Codon_ID': row['Codon_ID'],
                    'Position': pos,
                    'Alt': alt,
                    'Ref': ref,
                    'Reads': row['Reads']
                }
                # Append the new row to the new_rows list
                mutation_list.append(new_row)

    # Create a new DataFrame from the new_rows list
    df_new_mutations = pd.DataFrame(mutation_list) 
    df_new_mutations.to_csv('resmarker_new_mutations.txt', sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--allele_data_path", required=True, help="Path to allele_data.txt")
    parser.add_argument("--res_markers_info_path", required=True, help="Path to resistance marker table")
    parser.add_argument("--refseq_path", required=True, help="Path to reference sequences [fasta]")
    parser.add_argument("--n-cores", type=int, default=1, help="Number of cores to use")
    args = parser.parse_args()
    main(args)