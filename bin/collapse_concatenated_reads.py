#!/usr/bin/env python3
"""
Collapse concatenated reads.

Usage:
    python collapse_concatenated_reads.py --clusters path/to/asv_table.tsv
"""

import argparse
import pandas as pd
import sys

# TODO: add unit tests using nftest
# TODO: test full pipeline 
# TODO: remove old R script

def main():
    parser = argparse.ArgumentParser(description='Collapse concatenated reads with 10 Ns')
    parser.add_argument('--clusters', required=True, help='Path to clustered ASV table file')
    parser.add_argument('--output', help='Output file path (default: clusters.concatenated.collapsed.txt)')
    
    args = parser.parse_args()
    
    try:
        # Read the ASV table
        print(f"Reading clusters from: {args.clusters}", file=sys.stderr)
        df = pd.read_csv(args.clusters, sep='\t')
        
        # Check if required columns exist
        required_columns = ['asv', 'sampleID', 'locus', 'reads']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f"Error: Missing required columns: {missing_columns}", file=sys.stderr)
            print(f"Available columns: {list(df.columns)}", file=sys.stderr)
            sys.exit(1)
        
        # Filter for ASVs with exactly 10 Ns
        print("Processing ASVs with exactly 10 Ns...", file=sys.stderr)
        concatenated_df = df[df['asv'].str.count('N') == 10]
        unconcatenated_df = df[df['asv'].str.count('N') != 10]
        
        print(f"Original table: {len(df)} rows", file=sys.stderr)
        print(f"Concatenated ASVs: {len(concatenated_df)} rows", file=sys.stderr)
        print(f"Non-concatenated ASVs: {len(unconcatenated_df)} rows", file=sys.stderr)
        
        # Process concatenated ASVs if any exist
        if len(concatenated_df) > 0:
            collapsed_df = collapse_asvs(concatenated_df)
            print(f"Collapsed concatenated ASVs: {len(collapsed_df)} rows", file=sys.stderr)
    
            # Combine with non-concatenated ASVs (leave unconcatenated rows untouched)
            final_df = pd.concat([unconcatenated_df, collapsed_df], ignore_index=True)
            
            # Calculate normalized reads and allele counts without re-aggregating non-concatenated ASVs
            grp = final_df.groupby(['sampleID', 'locus'], sort=False)
            final_df = final_df.copy()
            final_df['norm.reads.locus'] = final_df['reads'] / grp['reads'].transform('sum')
            final_df['n.alleles'] = grp['asv'].transform('count')

        else:
            final_df = unconcatenated_df

        # Output results with all columns ordered correctly
        output_file = args.output if args.output else "clusters.concatenated.collapsed.txt"
        desired_cols = ['sampleID', 'locus', 'asv', 'reads', 'allele', 'norm.reads.locus', 'n.alleles']
        final_df = final_df[desired_cols]
        final_df.to_csv(output_file, sep='\t', index=False)
        print(f"Output saved to: {output_file}", file=sys.stderr)
        print(f"Final table: {len(final_df)} rows", file=sys.stderr)
            
    except FileNotFoundError:
        print(f"Error: File '{args.clusters}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

def collapse_asvs(concat_clusters):
    # Work on a copy to avoid SettingWithCopyWarning
    concat_clusters = concat_clusters.copy()
    concat_clusters.loc[:, 'left'] = concat_clusters['asv'].str.split('NNNNNNNNNN').str[0]
    concat_clusters.loc[:, 'right'] = concat_clusters['asv'].str.split('NNNNNNNNNN').str[1]
    concat_clusters.loc[:, 'len_left'] = concat_clusters['left'].str.len()
    concat_clusters.loc[:, 'len_right'] = concat_clusters['right'].str.len()
    # Avoid groupby.apply deprecation by iterating groups explicitly
    outputs = []
    for (_, _), g in concat_clusters.groupby(['sampleID', 'locus'], sort=False):
        outputs.append(collapse_group(g))
    if len(outputs) == 0:
        return pd.DataFrame(columns=['sampleID','locus','asv','reads','allele'])
    return pd.concat(outputs, ignore_index=True)
     
def collapse_group(df):
    """Collapse sequences within a group by finding common prefixes and suffixes."""
    options = []
    for _, row in df.iterrows():
        options.append({'left': row.left, 'right': row.right})

    new_lefts = []
    new_rights = []
    
    for _, row in df.iterrows(): 
        new_left = row.left
        new_right = row.right
        
        for option in options: 
            # Check if sequences can be collapsed based on prefix/suffix matching
            left_match = (row.left.startswith(option['left']) or 
                         option['left'].startswith(row.left))
            right_match = (row.right.endswith(option['right']) or 
                          option['right'].endswith(row.right))
            
            if left_match and right_match:
                # Take the shorter sequence for collapsing
                if len(option['left']) < len(new_left):
                    new_left = option['left']
                if len(option['right']) < len(new_right):
                    new_right = option['right']
                    
        new_lefts.append(new_left)
        new_rights.append(new_right)
    
    df['new_left'] = new_lefts
    df['new_right'] = new_rights
    
    # Group by the collapsed sequences and sum reads
    df_new = df.groupby(['sampleID', 'locus', 'new_left', 'new_right']).reads.sum().reset_index()
    df_new['asv'] = df_new.new_left + 'NNNNNNNNNN' + df_new.new_right
    df_new = df_new[['sampleID', 'locus', 'asv', 'reads']]
    
    # Add allele column
    allele_df = df[['locus','asv','allele']].drop_duplicates()
    df_new = df_new.merge(allele_df, on=['locus','asv'], how='left')
    # Convert concat to single N 
    df_new['asv'] = df_new.asv.str.replace('NNNNNNNNNN','N')
    return df_new

if __name__ == "__main__":
    main()