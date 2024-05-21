#!/usr/bin/env python3

"""
Created: May 16, 2024
Updated: May 16, 2024
Author(s): Todd Lenz, tlenz001@ucr.edu

Plots TSS coverage.
"""

import sys
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from  matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages

def parse_args(args):
    parser = argparse.ArgumentParser(description='Check the help flag')
    parser.add_argument('-b',
                        '--bed',
                        dest='bed',
                        help='REQUIRED: bed file(s) to plot.',
                        required=True)
    parser.add_argument('-g',
                        '--gff',
                        dest='gff',
                        help='Gene data in GFF format. Coding regions are '
                             'used to determine which coordinates within BED '
                             'files to plot.',
                        required=True,
                        default=None)
    parser.add_argument('-l',
                        '--gene_list',
                        dest='gene_list',
                        help='List of genes to plot, separated by spaces, '
                             'commas, newlines or tabs. Genes will be '
                             'extracted from the provided GFF file.',
                        required=False,
                        default=None)
    parser.add_argument('-s',
                        '--samples',
                        dest='samples',
                        help='Sample name(s)',
                        required=False)
    parser.add_argument('-o',
                        '--output',
                        dest='output',
                        help='Output directory.',
                        required=False)
    return parser.parse_args()

def input_params(args):
    """
    Parse input arguments and load bed files into dataframes. If 
    provided, genes will be extracted from an input gff file.
    """
    if args.bed:
        bed = pd.read_csv(args.bed, sep='\t', header=None)
        if not args.output:
            out = ''.join([os.path.dirname(os.path.abspath(args.bed)), 
                           '/gene_barplot.pdf'])
        else:
            out = args.output
    else:
        print("ERROR: No input BED file found. A BED file is required using the "
              "-b flag.")
        sys.exit()
    if args.gff:
        genes = filter_gff(pd.read_csv(args.gff, sep='\t', header=None, 
                                       usecols=[0,2,3,4,6,8]))
        if args.gene_list:
            genes = extract_genes(genes, args.gene_list)
        genes.columns = range(len(genes.columns))
    else:
        print("No GFF file found. A GFF file is required using the -g flag.")
        sys.exit()
    return bed, genes, out

def filter_gff(gff):
    """ Extracts gene accessions, names and descriptions. """
    genes = gff[(gff[2] == 'protein_coding_gene') | (gff[2] == 'ncRNA_gene') | 
                (gff[2] == 'pseudogene') | (gff[2] == 'gene')
                ].reset_index(drop=True).copy()
    genes[9] = genes[8].str.extract(r'ID=(.*?);', expand=True)
    genes[10] = genes[8].str.extract(r'Name=(.*?);', expand=True)
    genes[11] = genes[8].str.extract(r'description=(.*?);', expand=True)
    genes.drop([2,8], axis=1, inplace=True)
    genes[11] = genes[11].str.replace(r'%2C', ',', regex=True)
    return genes

def extract_genes(genes, gene_list):
    """
    Filters genes by list provided in text file. Gene names can be 
    delimited by spaces, commas, newlines or tabs.
    """
    delimiter = check_delimiter(gene_list)
    list_of_genes = []
    with open(gene_list, 'r') as f:
        for line in f:
            items = line.strip().split(delimiter)
            list_of_genes.extend(items)
    genes = genes[genes[9].isin(list_of_genes)]
    return genes

def check_delimiter(gene_list):
    """ Identify delimiter for list of genes. """
    with open(gene_list, 'r') as f:
        lines = [f.readline().strip() for _ in range(5)]
    space_count = sum(line.count(' ') for line in lines)
    comma_count = sum(line.count(',') for line in lines)
    newline_count = sum(line.count('\n') for line in lines)
    tab_count = sum(line.count('\t') for line in lines)
    max_count = max(space_count, comma_count, newline_count, tab_count)
    if space_count == max_count:
        return ' '
    if comma_count == max_count:
        return ','
    elif newline_count == max_count:
        return '\n'
    elif tab_count == max_count:
        return '\t'
    else:
        print('Can\'t determine delimiter of gene_list.')

def get_read_counts(bed, genes, fp, tp):
    read_counts = []
    for gene in genes.itertuples():
        if gene[4] == '+':
            bed_subset = bed[(bed[0] == gene[1]) & (bed[1] >= gene[2] - fp) & (bed[1] <= gene[2] + tp)]
            bed_subset.index = range(-fp, tp + 1)
            read_counts.append(bed_subset)
        else:
            bed_subset = bed[(bed[0] == gene[1]) & (bed[1] >= gene[3] - tp) & (bed[1] <= gene[3] + fp)]
            bed_subset.index = range(-fp, tp + 1)
            read_counts.append(bed_subset)
    return read_counts

def plot_read_counts(read_counts, fp, tp, out):
    pdf = PdfPages(out)
    plt.plot(read_counts)
    plt.xlabel('Position relative to TSS (bp)')
    plt.ylabel('Average Read Count')
    plt.title('Read Count Distribution Around TSS')
    pdf.savefig()
    plt.close()
    pdf.close()

def main():
    args = parse_args(sys.argv[1:])
    bed, genes, out = input_params(args)
    fp = 500
    tp = 1500
    read_counts = get_read_counts(bed, genes, fp, tp)
    merged_read_counts = pd.concat([df[[2]] for df in read_counts], axis=1)
    merged_read_counts = merged_read_counts.mean(axis=1)
    plot_read_counts(merged_read_counts, fp, tp, out)

if __name__ == '__main__':
    main()
