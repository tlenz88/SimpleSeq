#!/usr/bin/env python3

## Created: February 10, 2022
## Updated: May 17, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Converts a BED file to a WIG file for use with IGV.

import sys
import csv

def bed_to_wig(bed_file):
    if len(sys.argv) != 2:
        print("Usage: python3 bed2wig.py input.bed")
        sys.exit(1)
    
    bed_to_wig(bed_file)
    wig_file = bed_file.replace('.bed', '.wig')
    read_num = ['track type=wiggle_0 name=track_label']
    
    with open(bed_file, 'r') as bed:
        for row in csv.reader(bed, delimiter='\t'):
            if row[1] == '1':
                read_num.append(f'fixedStep chrom={row[0]} start=1 step=1')
            read_num.append(row[2])
    
    with open(wig_file, 'w') as wig:
        wig.write('\n'.join(read_num))

if __name__ == "__main__":
    bed_to_wig()