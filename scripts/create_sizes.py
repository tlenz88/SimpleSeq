#!/usr/bin/env python3

"""
Created: February 10, 2022
Updated: March 3, 2024
Author(s): Todd Lenz, tlenz001@ucr.edu

Finds the length of all chromosomes in a FASTA file and outputs a tab-delimited
text file with two columns: chromosome and length. The output will be stored in
the same directory as the input FASTA file.
"""

import sys
from Bio import SeqIO
import os

input_fasta = SeqIO.parse(open(sys.argv[1]), "fasta")
outfile = open("".join([os.path.dirname(os.path.abspath(sys.argv[1])), "/", 
                        os.path.splitext(os.path.basename(sys.argv[1]))[0], 
                        ".chrom.sizes"]),"w")
for chrom in input_fasta:
    outfile.write(chrom.id + "\t" + str(len(chrom)) + "\n")
outfile.close()
