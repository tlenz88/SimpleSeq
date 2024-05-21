#!/usr/bin/env python3

## Created: February 10, 2022
## Updated: May 17, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

## Converts a GAF file to a MAP file.

import sys
import csv
from collections import defaultdict

gene_to_GO = defaultdict(set)

with open(sys.argv[1], 'r') as file:
    csv_reader = csv.reader(file, delimiter='\t')
    for row in csv_reader:
        if len(row) >= 17:
            gene = row[1]
            GO = row[4]
            gene_to_GO[gene].add(GO)

with open(sys.argv[1].replace('.gaf', '.map'), 'w') as map:
    for gene, GO_set in gene_to_GO.items():
        GO_str = '\t'.join(sorted(GO_set))
        map.write(f"{gene}\t{GO_str}\n")