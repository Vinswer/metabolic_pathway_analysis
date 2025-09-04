#!/usr/bin/env python3
"""
Description: Merges gene counts
Usage: python3 <count_filepath> <output_filepath>
    count_filepath: filepath of input counts file, one or more
    output_filepath: filepath of output merged counts file
"""

from sys import argv
import pandas as pd

def merge_counts(count_filepaths, output_filepath):
    merged_df = pd.DataFrame(columns=['gene_id'])
    for count_file in count_filepaths:
        df = pd.read_csv(count_file)
        merged_df = pd.merge(merged_df, df, on='gene_id', how='outer')

    columns = ['gene_id'] + sorted([col for col in merged_df.columns if col != 'gene_id'])
    df_sorted = merged_df[columns]
    output_file = f"{output_filepath}"
    df_sorted.to_csv(output_file, index=False)

def main():
    """Main function"""
    count_filepaths = []
    for filepath in argv[1:-1]:
        count_filepaths.append(filepath)
    output_filepath = argv[-1]

    merge_counts(count_filepaths, output_filepath)

if __name__ == "__main__":
    main()
