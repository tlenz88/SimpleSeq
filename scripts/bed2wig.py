#!/usr/bin/env python3

# Used to convert a BED file to a WIG file for use with IGV

import sys
import csv

wig = (str(sys.argv[1])).replace(".bed", ".wig")
read_num = []
read_num.append("track type=wiggle_0 name=track_label")

for i in csv.reader(open(sys.argv[1], 'r'), delimiter='\t'):
	if i[1] == "1":
		read_num.append("fixedStep chrom=%s start=1 step=1" % i[0])
	read_num.append(i[2])

with open(wig, 'w+') as f:
	f.write("%s" % "\n".join(read_num))
