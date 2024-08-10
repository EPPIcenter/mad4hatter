#!/usr/bin/env python3
import argparse
import os
import pandas as pd

# TODO: Change to dictionary


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
    concatenated_df = concatenated_df.groupby(['amplicon', 'amplicon_start', 'amplicon_end', 'ampInsert_start', 'ampInsert_end', 'rev_primer',
                                               'amplicon_length', 'ampInsert_length', 'fwd_primer'])['pool'].agg(lambda x: ','.join(x)).reset_index()
    concatenated_df.sort_values('amplicon', inplace=True)
    concatenated_df.reset_index(inplace=True, drop=True)
    return concatenated_df


def build_amplicon_info():
    args = parse_args_build_amplicon_info()
    amplicon_info_df = concatenate_tables(args.amplicon_info_paths, args.pools)
    amplicon_info_df.to_csv(
        args.amplicon_info_output_path, index=False, sep='\t')


if __name__ == "__main__":
    build_amplicon_info()
