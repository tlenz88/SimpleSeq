#!/usr/bin/env python

import sys
import csv

lst = dict()

with open(sys.argv[1], "r") as f:
	csv = csv.reader(f, delimiter="\t")
	for i in csv:
		if len(i) >= 17:
			gene = i[1]
			GO = i[4]
			if gene not in lst.keys():
				lst[gene] = "\t" + GO
			else:
				if GO not in lst[gene]:
					lst[gene] += ", " + GO

with open(sys.argv[1].replace(".gaf", ".map"), "w+") as f:
	for key in lst.keys():
		f.write("%s%s\n" % (key, lst[key]))