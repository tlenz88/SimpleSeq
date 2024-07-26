#!/usr/bin/env python3

## Created: March 16, 2024
## Updated: May 17, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Performs count-per-million normalization of BED file.

import sys
import numpy as np
import pandas as pd
import os

def cpm_normalize():
    if len(sys.argv) != 2:
        print("Usage: python3 cpm_normalization.py input.bed")
        sys.exit(1)

    input_file = sys.argv[1]
    df = pd.read_csv(input_file, sep='\t', header=None)
    mmr = df[2].sum() / 1e6
    df[2] = df[2].div(mmr)
    base_name = os.path.basename(input_file)
    file_name, file_ext = os.path.splitext(base_name)
    output_file = os.path.join(os.path.dirname(input_file), 
                            f"{file_name}_cpm{file_ext}")
    df.to_csv(output_file, sep='\t', header=False, index=False)
    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    cpm_normalize()
