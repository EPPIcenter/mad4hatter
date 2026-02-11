#!/usr/bin/env python3
import argparse
import pandas as pd


def parse_args_build_resmarker_info():
    parser = argparse.ArgumentParser()
    parser.add_argument('--amplicon_info', type=str, required=True,
                        help='Path to file containing panel information')
    parser.add_argument('--principal_resmarkers', type=str, required=True,
                        help='Full list of all resmarkers of interest')
    parser.add_argument('--resmarker_info_output_path',
                        type=str, default='resmarker_info.tsv')
    return parser.parse_args()


def build_resmarker_info():
    args = parse_args_build_resmarker_info()

    panel_info_df = pd.read_csv(args.amplicon_info, dtype={
                                'GeneID': str}, sep='\t')
    full_resmarker_df = pd.read_csv(args.principal_resmarkers, dtype={
        'GeneID': str}, sep='\t')

    data = {'GeneID': [], 'Gene': [], 'CodonID': [], 'chrom': [], 'start': [], 'stop': [
    ], 'strand': [], 'target_name': [], 'CodonStart': []}

    for _, row in full_resmarker_df.iterrows():
        chr = row.chr
        codon_start = row.start
        codon_stop = row.stop

        panel_overlap = panel_info_df[(panel_info_df.chrom==chr) & (
            panel_info_df.insert_start <= codon_start) & (panel_info_df.insert_end >= codon_stop)]
        if len(panel_overlap) >= 1:
            for _, overlap_row in panel_overlap.iterrows():
                data['GeneID'].append(row.GeneID)
                data['Gene'].append(row.Gene)
                data['CodonID'].append(row.CodonID)
                data['chrom'].append(row.chr)
                data['start'].append(row.start)
                data['stop'].append(row.stop)
                data['strand'].append(row.strand)
                data['target_name'].append(overlap_row.target_name)
                codon_start = row.start-overlap_row.insert_start
                data['CodonStart'].append(codon_start)
    resmarker_df = pd.DataFrame(data=data)
    resmarker_df.drop_duplicates(inplace=True)
    resmarker_df.reset_index(inplace=True, drop=True)
    if not resmarker_df.empty:
        resmarker_df.to_csv(
            args.resmarker_info_output_path, index=False, sep='\t')


if __name__ == "__main__":
    build_resmarker_info()
