#!/usr/bin/env python
"""
Created June 13, 2022
Updated August 8, 2023

@author: Todd Lenz, tlenz001@ucr.edu

Subtracts control reads (input/IGG) from conditional reads.

Required input arguments:
1. b1--BED file containing genome-wide per-base coverage for 
   conditional sample. 
   To generage this file from a BAM file:
   bedtools genomecov -d -ibam ${i}_sorted.bam > ${i}_sorted.bed
2. b2--BED file containing genome-wide per-base coverage for 
   control sample.

Optional input arguments:
1. normalize--To normalize data in BED files so that samples can be 
   directly compared, indicate a method for normalization. Data can be 
   normalized using counts-per-million (CPM). More methods will be 
   added in the future.
2. outdir--Path to output directory. File name will be extracted from 
   conditional input. If no output directory is given, the output will 
   be the saved in the same directory as the conditional input.
"""


import sys
import argparse
import pandas as pd
import os


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b1',
                        dest='bed1',
                        help='BED file of target sample.',
                        required=True)
    parser.add_argument('-b2',
                        dest='bed2',
                        help='BED file of input or IGG.',
                        required=True)
    parser.add_argument('-n',
                        '--normalize',
                        dest='normalize',
                        help='Normalize read counts for each sample to '
                             'account for differences in sequencing depth. '
                             'Data can be normalized using CPM (counts-per-'
                             'million) normalization.',
                        required=False,
                        choices=['CPM'],
                        default=None)
    parser.add_argument('-o',
                        '--output',
                        dest='output',
                        help='Path to output directory. If no output '
                             'directory is given, the output will be saved '
                             'in the same directory as the conditional '
                             'input.',
                        required=False)
    return parser.parse_args()


def normalize_df(df, norm):
    """ Normalize data from BED files. """
    if norm == 'CPM':
        for i in range(2,len(df.columns)):
            mmr = df[i].sum() / 1000000
            df[i] = df[i].div(mmr)
    else:
        print('Error: Provided method of normalization is not valid.')
        pass
    return df


def subtract_chip(bed1, bed2):
    """ Subtract control reads from target reads. """
    bed = bed1.merge(bed2, on=[0,1]).fillna(0)
    bed.columns = range(len(bed.columns))
    bed[4] = bed[3] - bed[2]
    bed.drop([2,3], axis=1, inplace=True)
    bed.loc[bed[4] < 0, 4] = 0
    bed.columns = range(len(bed.columns))
    return bed


def main():
    args = parse_args(sys.argv[1:])
    bed1 = pd.read_csv(args.bed1, sep='\t', header=None)
    bed2 = pd.read_csv(args.bed2, sep='\t', header=None)
    if args.normalize:
        print('Normalizing data using %s' % str(args.normalize).upper())
        norm_bed1 = normalize_df(bed1, args.normalize)
        norm_bed2 = normalize_df(bed2, args.normalize)
        diff_bed = subtract_chip(norm_bed1, norm_bed2)
    else:
        diff_bed = subtract_chip(bed1, bed2)
    if args.output:
        outfile = ''.join([os.path.splitext(os.path.basename(args.bed1))[0], 
                           '_diff.bed'])
        diff_bed.to_csv('%s/%s' % (args.output, outfile), sep = '\t', 
                        header = False, index = False)
    else:
        outdir = os.path.dirname(os.path.abspath(args.bed1))
        outfile = ''.join([os.path.splitext(os.path.basename(args.bed1))[0], 
                           '_diff.bed'])
        diff_bed.to_csv('%s/%s' % (outdir, outfile), sep = '\t', 
                        header = False, index = False)


if __name__ == '__main__':
    main()
