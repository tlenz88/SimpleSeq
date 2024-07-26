#!/usr/bin/env python3

"""
Created: July 10, 2023
Updated: January 30, 2024
Author(s): Todd Lenz, tlenz001@ucr.edu

Generates coverage tracks for all genes in a gff file or those in a
provided list. Track lengths are normalized with respect to the longest
plotted track--i.e. the longest gene will fill the width of the figure
and all other genes are proportionally plotted against it.
"""

import sys
import os
import math
import re
import argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-b',
                        '--bed',
                        dest='bed',
                        help='Tab-delimited BED files with per-base '
                             'read coverage (bedtools genomecov -d).',
                        nargs='+',
                        required=True)
    parser.add_argument('-g',
                        '--gff',
                        dest='gff',
                        help='Gene data in GFF format. Coding regions are '
                             'used to determine which coordinates within BED '
                             'files to plot and exons are used to generate '
                             'arrows to annotate the x-axis of each plot.',
                        required=True)
    parser.add_argument('-s',
                        '--samples',
                        dest='samples',
                        help='Sample name(s) used to annotate barplot(s).',
                        nargs='+',
                        required=False)
    parser.add_argument('-o',
                        '--output',
                        dest='output',
                        help='Output file path and name. If no file is given, '
                             'the output will be saved in the same directory '
                             'as the first BED file with a default file name.',
                        required=False)
    parser.add_argument('-r',
                        '--resolution',
                        dest='resolution',
                        help='Resolution at which to bin read count data. '
                             'Binning will begin at coordinate 1 and continue '
                             'to the end of the gene, therefore the last bin '
                             'may be less than the given value.',
                        required=False,
                        default=10,
                        type=int)
    parser.add_argument('-l',
                        '--gene_list',
                        dest='gene_list',
                        help='Text file containing a list of genes to plot. '
                             'Genes can be separated by spaces, commas, '
                             'newlines or tabs. Genes will be extracted from '
                             'the provided GFF file.',
                        nargs='+',
                        required=False)
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
    parser.add_argument('-d',
                        '--distance',
                        dest='distance',
                        help='Distance (bp) from the ends of genes--5\' from '
                             'the TSS and 3\' from the end of the gene--that '
                             'will be plotted in the figure.',
                        required=False,
                        default=500,
                        type=int)
    parser.add_argument('--global_ymax',
                        dest='global_ymax',
                        action='store_true',
                        help='Use the same maximum y-value for all plots. If not set, each plot will have its own y-max.')
    return parser.parse_args()

def input_params(args):
    """
    Parse input arguments, loading gff and bed files into dataframes 
    and extract genes and exons from gff. If provided, the samples, 
    resolution, and distance arguments will also be parsed.
    """
    genes, exons = filter_gff(pd.read_csv(args.gff, sep='\t', header=None, 
                                          usecols=[0,2,3,4,6,8]))
    if args.gene_list:
        if len(args.gene_list) == 1:
            for gl in args.gene_list:
                if os.path.exists(gl):
                    genes, exons = extract_genes(genes, exons, gl)
                else:
                    genes = genes[genes[9] == gl]
        else:
            genes = genes[genes[9].isin(args.gene_list)]
    samples = []
    for i in args.bed:
        bed = pd.read_csv(i, sep='\t', header=None)
        try:
            df = df.merge(bed, on=[0,1])
        except:
            df = bed
            if not args.output:
                out = ''.join([os.path.dirname(i), '/ChIP_barplot.pdf'])
            else:
                out = args.output
        if not args.samples or len(args.samples) != len(args.bed):
            samples.append(os.path.basename(i)[:-4])
    df.columns = range(len(df.columns))
    if len(samples) == 0:
        samples = args.samples
    if args.normalize:
        df = normalize_df(df, args.normalize)
    res = int(args.resolution)
    dist = int(args.distance)
    return genes, exons, df, out, samples, res, dist, args.global_ymax

def filter_gff(gff):
    """ Extracts gene and exon accessions, names and descriptions. """
    exons = gff[gff[2] == 'exon'].drop([2], axis=1).copy()
    genes = gff[(gff[2] == 'protein_coding_gene') | (gff[2] == 'ncRNA_gene') | 
                (gff[2] == 'pseudogene') | (gff[2] == 'gene')
                ].reset_index(drop=True).copy()
    genes[9] = genes[8].str.extract(r'ID=(.*?);', expand=True)
    genes[10] = genes[8].str.extract(r'Name=(.*?);', expand=True)
    genes[11] = genes[8].str.extract(r'description=(.*?);', expand=True)
    genes.drop([2,8], axis=1, inplace=True)
    #genes[11] = genes[11].str.replace(r'%2C', ',', regex=True)
    exons[9] = exons[8].str.extract(r'exon_(.*?)\.[1-9]-E', expand=True)
    exons[10] = exons[8].str.extract(r'exon_.*?\.([1-9]).*?;', expand=True)
    exons = exons[exons[10] == '1'].reset_index(drop=True).drop([8,10],axis=1)
    return genes, exons

def extract_genes(genes, exons, gene_list):
    """
    Filters genes and exons by list of genes provided in text file. 
    Gene names can be delimited by spaces, commas, newlines or tabs.
    """
    delimiter = check_delimiter(gene_list)
    list_of_genes = []
    with open(gene_list, 'r') as f:
        for line in f:
            items = line.strip().split(delimiter)
            list_of_genes.extend(items)
    genes = genes[genes[9].isin(list_of_genes)]
    exons = exons[exons[9].isin(list_of_genes)]
    return genes, exons

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
        print("Can't determine delimiter of gene_list.")

def normalize_df(df, norm):
    """ Normalize data from BED files. """
    if norm == 'CPM':
        for i in range(2,len(df.columns)):
            mmr = df[i].sum() / 1000000
            df[i] = df[i].div(mmr)
    else:
        print("Provided method of normalization is not valid so data will not "
              "be normalized. Run script with -h flag to see valid "
              "normalization methods.")
        pass
    return df

def extract_gene_regions(df, genes, dist):
    """ Extract read counts for genes of interest. """
    for g in genes.itertuples():
        gene_df = df[(df[0] == g[1]) & (df[1] >= g[2] - dist) & 
                     (df[1] <= g[3] + dist)].reset_index(drop=True).copy()
        gene_df[0] = g[5]
        gene_df[1] = [*range(1, len(gene_df) + 1)]
        if pd.isna(g[6]):
            gene_df[len(gene_df.columns)] = ''.join(
                [g[5], ' - ', g[7]])
        else:
            gene_df[len(gene_df.columns)] = ''.join(
                [g[5], ' - ', g[6], ': ', g[7]])
        gene_df[len(gene_df.columns)] = g[4]
        try:
            gdf = pd.concat([gdf, gene_df], ignore_index=True)
        except:
            gdf = gene_df
    return gdf

def data_binning(df, res):
    """ Bins read counts per gene. """
    df[len(df.columns)] = (df.iloc[:, 1] - 1) // res * res + 1
    df = df[list(df.columns)[0:2] + 
            list(df.columns)[-3:] + 
            list(df.columns)[2:-3]]
    df.columns = range(len(df.columns))
    binned_df = df.groupby([0,2,3,4], 
                           as_index=False)[list(df.columns)[4:]].sum()
    binned_df.columns = range(len(binned_df.columns))
    return binned_df

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

def create_arrow(ax, exons, strand, res, max_yval, dist):
    """ Annotate each plot with arrows representing exons. """
    arrows = []
    start = dist // res + 1
    prev_end = 0
    prev_bin = 0
    for e in exons.itertuples():
        if e[2] % res != 0:
            start_bin = e[2] // res + 1
        else:
            start_bin = e[2] // res
        if e[3] % res != 0:
            end_bin = e[3] // res + 1
        else:
            end_bin = e[3] // res
        if e[0] == 0 or len(arrows) == 0:
            arrows.append([start, start + end_bin - start_bin])
            prev_bin = end_bin - start_bin + start
            prev_end = end_bin
        else:
            arrows.append([start_bin - prev_end + prev_bin, 
                           end_bin - prev_end + prev_bin])
            prev_bin = end_bin - prev_end + prev_bin
            prev_end = end_bin
    if strand == '+':
        for arrow in arrows:
            ax.arrow(arrow[0], -max_yval * .5, arrow[1] - arrow[0], 0, 
                     width=max_yval/4, head_width=max_yval/2, head_length=10, 
                     facecolor='#FF0000', edgecolor='#FF0000', 
                     length_includes_head=True, clip_on=False)
    elif strand == '-':
        for arrow in arrows[::-1]:
            ax.arrow(arrow[1], -max_yval*.5, -(arrow[1] - arrow[0]), 0, 
                     width=max_yval/4, head_width=max_yval/2, head_length=10, 
                     facecolor='#FF0000', edgecolor='#FF0000', 
                     length_includes_head=True, clip_on=False)
    return ax

def main():
    args = parse_args(sys.argv[1:])
    genes, exons, df, out, samples, res, dist, global_ymax = input_params(args)
    df = extract_gene_regions(df, genes, dist)
    df = data_binning(df, res)
    max_yval = df[list(df.columns)[4:]].max().max() if global_ymax else None
    max_xval = max(df[0].value_counts())
    pdf = PdfPages(out)
    for gene, idx in zip(df[0].unique(), range(len(df[0].unique()))):
        gene_df = df[df[0] == gene].copy()
        gene_df[3] = [x + 1 for x in list(range(len(gene_df)))]
        fig = plt.figure()
        fig.set_figheight(len(samples))
        fig.set_figwidth(20)
        sample_colors = ['#F3766E', '#1CBDC2']
        local_max_yval = gene_df[list(gene_df.columns)[4:]].max().max() if not global_ymax else max_yval
        for i in [*range(len(samples))]:
            ax = plt.subplot2grid((len(samples)+1, 20), (i, 0), 
                                  colspan=int(math.ceil(max(gene_df[3]) / 
                                                        max_xval * 20)), 
                                  rowspan=1)
            barplt = plt.bar(np.array(gene_df[3]), np.array(gene_df[4+i]), 
                             width=1, color=sample_colors[i])
            ax = custom_ax_params(ax, samples[i], local_max_yval, max(gene_df[3]))
            if i == len(samples) - 1:
                ax = create_arrow(ax, exons[exons[9] == gene_df[0].iloc[0]
                                            ].reset_index(drop=True), 
                                            gene_df[2].iloc[0], res, 
                                            local_max_yval, dist)
        fig.suptitle(gene_df[1].iloc[0])
        fig.subplots_adjust(hspace=0)
        fig.tight_layout(h_pad=0)
        pdf.savefig()
        plt.close()
    pdf.close()

if __name__ == '__main__':
    main()
