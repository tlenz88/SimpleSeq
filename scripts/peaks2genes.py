#!/usr/bin/env python3

import re
import os
import sys
import csv
import pandas as pd
import numpy as np


def assign_args(args):
    for i in args:
        if os.path.splitext(i)[1] == '.gff':
            gff = pd.read_csv(i, sep='\t', header=None)
            genes = gff.loc[(gff[2] == 'protein_coding_gene') | (gff[2] == 'ncRNA_gene')]
            accession = genes[8].str.extract(r'ID=(.*?);', expand=True)
            name = genes[8].str.extract(r'Name=(.*?);', expand=True)
            description = genes[8].str.extract(r'description=(.*?);', expand=True)
            genes = pd.concat([genes, accession, name, description], axis=1)
            genes.columns = range(genes.columns.size)
            genes.reset_index(drop=True, inplace=True)
            genes[11] = genes[11].str.replace(r'%2C', ',', regex=True)
            genes.sort_values([0,3], ascending=True, inplace=True)
            genes.reset_index(drop=True, inplace=True)
        elif os.path.splitext(i)[1] == '.narrowPeak' or os.path.splitext(i)[1] == '.broadPeak':
            peaks = pd.read_csv(i, sep='\t', header=None)
            peaks.sort_values([0, 1], ascending=True, inplace=True)
            peaks.reset_index(drop=True, inplace=True)
            outdir = os.path.dirname(os.path.abspath(i))
            outfile = ''.join([os.path.splitext(os.path.basename(i))[0], '_mapped.txt'])
        elif os.path.splitext(i)[1] == '.bed' or os.path.splitext(i)[1] == '.txt':
            peaks = pd.read_csv(i, sep='\t')
            peaks.sort_values(['chr', 'start'], ascending=True, inplace=True)
            peaks.reset_index(drop=True, inplace=True)
            outdir = os.path.dirname(os.path.abspath(i))
            outfile = ''.join([os.path.splitext(os.path.basename(i))[0], '_mapped.txt'])
        elif os.path.splitext(i)[1] == '.xlsx':
            peaks = pd.read_excel(i)
            peaks.sort_values(['chr', 'start'], ascending=True, inplace=True)
            peaks.reset_index(drop=True, inplace=True)
            outdir = os.path.dirname(os.path.abspath(i))
            outfile = ''.join([os.path.splitext(os.path.basename(i))[0], '_mapped.txt'])
        else:
            print()
            print("Input file %s is not a recognized input." % i)
            print("Run script with standard GFF and narrowPeak/broadPeak files.")
            print()
            sys.exit(0)
    return genes, peaks, outdir, outfile


def check_overlaps(sp, sg):
    overlap = 0
    sp_vals = list(range(sp[2],sp[3]))
    peak_overlap = []
    for g in sg.itertuples():
        g_vals = list(range(g[4],g[5]))
        ol = len(list(set(sp_vals) & set(g_vals)))
        if ol > overlap:
            overlap = ol
            peak_overlap = g[10:]
        else:
            continue
    return peak_overlap


def main():
    genes, peaks, outdir, outfile = assign_args(sys.argv[1:])
    gene_peaks = []
    for chrom in peaks['chr'].unique():
        sub_peaks = peaks[peaks['chr'] == chrom].copy()
        sub_genes = genes[genes[0] == chrom].copy()
        for sp in sub_peaks.itertuples():
            sg = sub_genes[(sub_genes[3] - 1000 < (sp[3] - sp[2]) / 2 + sp[2]) & (sub_genes[4] > (sp[3] - sp[2]) / 2 + sp[2])].copy()
            if len(sg.index) == 0:
                gene_peaks.append(['intergenic', '', ''])
            elif len(sg.index) == 1:
                gene_peaks.append(sg.iloc[:, 9:].values.flatten().tolist())
            else:
                gene_peaks.append(sg.iloc[:, 9:].values.flatten().tolist())
    gene_peaks_df = pd.DataFrame(gene_peaks)
    peaks_mapped = pd.concat([peaks, gene_peaks_df], axis=1)
    peaks_mapped.to_csv(''.join([outdir,'/',outfile]), sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()
