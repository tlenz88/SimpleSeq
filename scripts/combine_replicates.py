#!/usr/bin/env python3

## Created: March 16, 2024
## Updated: May 17, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Merges two BED files.

import sys
import pandas as pd

def combine_bed_replicates():
    if len(sys.argv) < 4:
        print("Usage: combine_replicates.py input1.bed "
              "input2.bed output.bed")
        sys.exit(1)

    try:
        df1 = pd.read_csv(sys.argv[1], sep='\t', header=None, 
                          names=['chrom', 'coord', 'reads_df1'])
        df2 = pd.read_csv(sys.argv[2], sep='\t', header=None, 
                          names=['chrom', 'coord', 'reads_df2'])
    except FileNotFoundError:
        print("One or more input files not found.")
        sys.exit(1)

    df = df1.merge(df2, on=['chrom', 'coord'])
    df['reads'] = df['reads_df1'] + df['reads_df2']
    df = df[['chrom', 'coord', 'reads']]
    df.to_csv(sys.argv[3], sep='\t', header=False, index=False)

if __name__ == "__main__":
    combine_bed_replicates()
