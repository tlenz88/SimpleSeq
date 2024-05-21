#!/usr/bin/env python3

"""
Created: June 13, 2022
Updated: December 12, 2023
Author(s): Todd Lenz, tlenz001@ucr.edu

Plots interaction heatmaps for Hi-C data.
"""

import sys
import os
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import font_manager
from matplotlib.backends.backend_pdf import PdfPages


def parse_args(args):
    parser = argparse.ArgumentParser(description='Check the help flag')
    parser.add_argument('-r',
                        '--read_counts',
                        dest='read_counts',
                        help='REQUIRED: tab-delimited matrix of reads per gene.',
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
                        '--sample',
                        dest='sample',
                        help='Sample name',
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
    if args.read_counts:
        read_counts = pd.read_csv(args.read_counts, sep='\t', index_col=0)
        read_counts.index.name = None
        if not args.output:
            out = ''.join([os.path.dirname(os.path.abspath(args.read_counts)), 
                           '/gene_barplot.pdf'])
        else:
            out = args.output
    else:
        print("No read count matrix found. A tab-delimited matrix is required "
              "using the -r flag.")
        sys.exit()
    if args.gff:
        genes = filter_gff(pd.read_csv(args.gff, sep='\t', header=None, 
                                       usecols=[0,2,3,4,8]))
        if args.gene_list:
            genes = extract_genes(genes, args.gene_list)
        genes.columns = range(len(genes.columns))
    else:
        print("No GFF file found. A GFF file is required using the -g flag.")
        sys.exit()
    return read_counts, genes, out


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


def main():
    args = parse_args(sys.argv[1:])
    read_counts, genes, out = input_params(args)

    mean_counts_sample1 = [np.mean(read_counts.loc[g, ['A3A', 'A3B']]) for g in genes[3]]
    mean_counts_sample2 = [np.mean(read_counts.loc[g, ['D2A', 'D2B']]) for g in genes[3]]

    err_sample1 = [np.min(read_counts.loc[g, ['A3A', 'A3B']]) for g in genes[3]]
    err_sample1 = [x - y for x, y in zip(mean_counts_sample1, err_sample1)]
    err_sample2 = [np.min(read_counts.loc[g, ['D2A', 'D2B']]) for g in genes[3]]
    err_sample2 = [x - y for x, y in zip(mean_counts_sample2, err_sample2)]

    pdf = PdfPages(out)
    fig, ax = plt.subplots(figsize=(12, 7))

    bar_width = 0.75
    gap_between_groups = bar_width
    gap_between_bars = 0.01
    x = np.arange(0, len(genes) * 2, 2)
    x1 = x - (bar_width + gap_between_bars) / 2
    x2 = x + (bar_width + gap_between_bars) / 2

    bars1 = plt.bar(x1, mean_counts_sample1, width=bar_width, 
                    color='#F3766E', label='WT')
    bars2 = plt.bar(x2, mean_counts_sample2, width=bar_width, 
                    color='#1CBDC2', label=r'$\Delta$V2')

    for bar, bar_error in zip(bars1, err_sample1):
        ax.errorbar(bar.get_x() + bar.get_width() / 2, bar.get_height(), 
                    yerr=[[bar_error], [bar_error]], color='black', 
                    capsize=bar_width * 1.5, capthick=bar_width / 1.5, 
                    elinewidth=bar_width / 1.25)

    for bar, bar_error in zip(bars2, err_sample2):
        ax.errorbar(bar.get_x() + bar.get_width() / 2, bar.get_height(), 
                    yerr=[[bar_error], [bar_error]], color='black', 
                    capsize=bar_width * 1.5, capthick=bar_width / 1.5, 
                    elinewidth=bar_width / 1.25)

    # Add some text for labels, title and custom x-axis tick labels, etc.
    font = {'family':'Times New Roman',
            'color':'black',
            'size':8}
    
    plt.ylabel('Normalized read count', labelpad=8, fontdict=font, ha='right')
    ax.set_ylim(bottom=0, top=5000)
    yticks = list(range(0, 5001, 1000))
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks, fontdict=font, ha='right')

    ax.spines['left'].set_linewidth(1)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['bottom'].set_position('zero')

    plt.xlabel('var genes', fontdict=font, ha='center')
    ax.set_xlim(left=-bar_width * 2, right=len(genes[3]) * 2)
    ax.set_xticks([i for i in range(0, len(genes[3]) * 2, 2)])
    ax.set_xticklabels(genes[3], fontdict=font, ha='center')
    plt.xticks(rotation=90, fontsize=8)

    ax.tick_params(axis='both', which='both', pad=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend(title='Sample', loc='upper center', 
               edgecolor='black', labelspacing=0.5, 
               fontsize=font['size'], title_fontsize=font['size'], 
               borderpad=0.5, prop={'family': font['family']}, 
               handlelength=1.5, handleheight=2)

    plt.tight_layout()
    pdf.savefig()
    plt.close()
    pdf.close()


if __name__ == '__main__':
    main()
