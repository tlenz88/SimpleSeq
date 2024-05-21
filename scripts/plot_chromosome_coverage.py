#!/usr/bin/env python3

"""
Created: June 13, 2022
Updated: December 12, 2023
Author(s): Todd Lenz, tlenz001@ucr.edu

Plots genomewide coverage for ChIP-seq data. Track lengths are
normalized with respect to the longest chromosome--i.e. the longest
chromosome will fill the width of the figure and all other chromosomes
are proportionally plotted against it.
"""


import sys
import os
import math
import re
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b',
                        '--bed',
                        dest='bed',
                        help='Tab-delimited BED files with per-base '
                             'read coverage.',
                        nargs='+',
                        required=True)
    parser.add_argument('-g',
                        '--gff',
                        dest='gff',
                        help='Gene data in GFF format. Coding regions are '
                             'used to determine which coordinates within BED '
                             'files to plot.',
                        required=False,
                        default=None)
    parser.add_argument('-s',
                        '--samples',
                        dest='samples',
                        help='Sample name(s) used to annotate barplot(s).',
                        nargs='+',
                        required=False,
                        default=None)
    parser.add_argument('-o',
                        '--output',
                        dest='output',
                        help='Output file path and name. If no file is '
                             'given, the output will be saved in the same '
                             'directory as the first BED file with a default '
                             'file name.',
                        required=False,
                        default=None)
    parser.add_argument('-r',
                        '--resolution',
                        dest='resolution',
                        help='Resolution at which to bin read count data. '
                             'Binning will begin at coordinate 1 and '
                             'continue to the end of the chromosome, '
                             'therefore the last bin may be less than the '
                             'given value.',
                        required=False,
                        default=10000,
                        type=int)
    parser.add_argument('-m',
                        '--centromeres',
                        dest='centromeres',
                        help='List of centromere coordinates.',
                        required=False,
                        default=None)
    parser.add_argument('-l',
                        '--gene_list',
                        dest='gene_list',
                        help='List of genes to plot, separated by spaces, '
                             'commas, newlines or tabs. Genes will be '
                             'extracted from the provided GFF file.',
                        required=False,
                        default=None)
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
    parser.add_argument('-y',
                        '--ymax',
                        dest='ymax',
                        help='Set the maximum y-value for all chromosomes.',
                        required=False,
                        default=None,
                        type=int)
    parser.add_argument('-c',
                        '--colors',
                        dest='colors',
                        help='List of colors to use for each sample.',
                        nargs='+',
                        required=False,
                        default=None)
    return parser.parse_args()


def input_params(args):
    """
    Parse input arguments and load bed files into dataframes. If 
    provided, genes will be extracted from an input gff file.
    """
    genes, centromeres = None, None
    if args.gff:
        genes = filter_gff(pd.read_csv(args.gff, sep='\t', header=None, 
                                       usecols=[0,2,3,4,8]))
        if args.gene_list:
            genes = extract_genes(genes, args.gene_list)
        genes.columns = range(len(genes.columns))
    if args.centromeres:
        centromeres = pd.read_csv(args.centromeres, sep='\t', header=None)
    samples = []
    for i in args.bed:
        bed = pd.read_csv(i, sep='\t', header=None)
        try:
            df = df.merge(bed, on=[0,1])
            df.columns = range(len(df.columns))
        except:
            df = bed
            if not args.output:
                out = ''.join([os.path.dirname(os.path.abspath(i)), '/ChIP_barplot.pdf'])
            else:
                out = args.output
        if not args.samples or len(args.samples) != len(args.bed):
            samples.append(os.path.basename(i)[:-4])
    if len(samples) == 0:
        samples = args.samples
    if args.normalize:
        df = normalize_df(df, args.normalize)
    if args.colors:
        sample_colors = args.colors
    else:
        sample_colors = sns.color_palette('colorblind', n_colors=len(args.bed))
    res = int(args.resolution)
    return genes, centromeres, df, out, samples, res, sample_colors


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


def normalize_df(df, norm):
    """ Normalize data from BED files. """
    if norm == 'CPM':
        for i in range(2,len(df.columns)):
            mmr = df[i].sum() / 1000000
            df[i] = df[i].div(mmr)
    else:
        print('Error: Provided method of normalization is not valid.')
        exit()
    return df


def data_binning(df, res):
    """ Bins read counts for each chromosome. """
    df[len(df.columns)] = df.iloc[:, 1] // res + 1
    df = df[list(df.columns)[0:1] + 
            list(df.columns)[-1:] + 
            list(df.columns)[2:-1]]
    df.columns = range(len(df.columns))
    df = df.groupby([0,1], as_index=False)[list(df.columns)[2:]].sum()
    return df


def custom_ax_params(ax, sample, max_yval, max_xval):
    """ Set font parameters and axes labels, limits and tick marks. """
    font = {'family':'serif',
            'color':'black',
            'weight':'bold',
            'size':10}
    ax.set_ylim(bottom=0, top=max_yval)
    plt.ylabel(sample, rotation='horizontal', labelpad=20, 
               fontdict=font, ha='right')
    ax.set_xlim(left=0, right=max_xval)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1)
    ax.spines['left'].set_capstyle('round')
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['bottom'].set_capstyle('round')
    ax.spines['bottom'].set_position('zero')
    return ax


def annotate_genes(ax, genes, res, max_yval):
    """ Annotate each plot with bars representing genes. """
    for g in genes.itertuples():
        ax.arrow(g[2] // res + 1, -max_yval * .2, (g[3]-g[2]) // res + 1, 
                 0, width=max_yval/4, head_width=0, head_length=0, 
                 facecolor='#FF0000', edgecolor='#FF0000', 
                 length_includes_head=True, clip_on=False)
    return ax


def annotate_centromeres(ax, centromeres, res, max_yval):
    """ Annotate each plot with bars representing genes. """
    for c in centromeres.itertuples():
        ax.arrow(c[2] // res + 1, -max_yval * .2, (c[3]-c[2]) // res + 1, 
                 0, width=max_yval/4, head_width=0, head_length=0, 
                 facecolor='#808080', edgecolor='#808080', 
                 length_includes_head=True, clip_on=False)
    return ax


def main():
    args = parse_args(sys.argv[1:])
    genes, centromeres, df, out, samples, res, sample_colors = input_params(args)
    df = data_binning(df, res)
    if args.ymax:
        #max_yval = df[list(df.columns)[2:]].max().max()
        max_yval = 5000
    max_xval = max(df[0].value_counts())
    pdf = PdfPages(out)
    num_plots = len(df[0].unique())*(len(samples)+1)
    plot_range = [*range(num_plots)]
    plot_idx = 0
    fig = plt.figure()
    fig.set_figheight(num_plots)
    fig.set_figwidth(20)
    H3K9_yval_list = [1370, 1968, 1912, 1906, 847, 1407, 1962, 2110, 1694, 1678, 1785, 1742, 1703, 666]
    #input_yval_list = [97, 99, 180, 151, 96, 99, 98, 103, 100, 115, 124, 192, 101, 98]
    y_num = 0
    for chr, idx in zip(df[0].unique(), range(len(df[0].unique()))):
        chr_df = df[df[0] == chr].copy()
        if args.ymax is None:
            #max_yval = chr_df[list(chr_df.columns)[2:]].max().max()
            max_yval = H3K9_yval_list[y_num]
            #max_yval = input_yval_list[y_num]
            y_num += 1
        for i in [*range(len(samples))]:
            ax = plt.subplot2grid((num_plots, 20), (plot_range[plot_idx], 0), 
                                  colspan=int(math.ceil(max(chr_df[1]) 
                                                        / max_xval * 20)), 
                                                        rowspan=1)
            barplt = plt.bar(np.array(chr_df[1]), np.array(chr_df[2+i]), 
                             width=1, color=sample_colors[i])
            ax = custom_ax_params(ax, samples[i], max_yval, max(chr_df[1]))
            plot_idx += 1
            if i == len(samples) - 1:
                plot_idx += 1
                if genes is not None:
                    ax = annotate_genes(ax, 
                                        genes[genes[0] == chr_df[0].iloc[0]
                                              ].reset_index(drop=True), 
                                              res, max_yval)
                if centromeres is not None:
                    ax = annotate_centromeres(ax, 
                                        centromeres[centromeres[0] == 
                                                    chr_df[0].iloc[0]
                                              ].reset_index(drop=True), 
                                              res, max_yval)
    fig.subplots_adjust(hspace=0)
    fig.tight_layout(h_pad=0)
    pdf.savefig()
    plt.close()
    pdf.close()


if __name__ == '__main__':
    main()
