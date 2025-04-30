#!/usr/bin/env python3
import argparse
import os
import pandas as pd


def parse_args_build_amplicon_info():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pools', type=str, required=True, nargs='+',
                        help='Pools corresponding to amplicon_info_paths')
    parser.add_argument('--amplicon_info_paths', type=str, required=True, nargs='+',
                        help='Paths to amplicon info tables for pools')
    parser.add_argument('--amplicon_info_output_path',
                        type=str, default='amplicon_info.tsv')
    return parser.parse_args()

def concatenate_tables(paths, pools):
    df_list = []
    for i in range(len(paths)):
        file = paths[i]
        pool = pools[i]
        if file == 'null':
            pass
        elif not os.path.exists(file):
            raise FileNotFoundError(file)
        else:
            df = pd.read_csv(file, sep='\t')
            df['pool'] = pool
            df_list.append(df)
    concatenated_df = pd.concat(df_list)
    concatenated_df = concatenated_df.groupby(['target_id', 'chrom', 'insert_start', 'insert_end', 'rev_primer','fwd_primer'])['pool'].agg(lambda x: ','.join(x)).reset_index()
    concatenated_df.sort_values('target_id', inplace=True)
    concatenated_df.reset_index(inplace=True, drop=True)
    concatenated_df['insert_length'] = concatenated_df.insert_end-(concatenated_df.insert_start+1) #TODO: change this to be calculated correctly 
    return concatenated_df


def build_amplicon_info():
    args = parse_args_build_amplicon_info()
    amplicon_info_df = concatenate_tables(args.amplicon_info_paths, args.pools)
    amplicon_info_df.to_csv(
        args.amplicon_info_output_path, index=False, sep='\t')


if __name__ == "__main__":
    build_amplicon_info()
