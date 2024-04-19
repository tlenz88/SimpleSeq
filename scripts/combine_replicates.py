#!/usr/bin/env python3

import sys
import pandas as pd

df1 = pd.read_csv(sys.argv[1], sep="\t", header=None)
df2 = pd.read_csv(sys.argv[2], sep="\t", header=None)
df = df1.merge(df2, on=[0,1])
df[2] = df['2_x'] + df['2_y']
df = df[[0,1,2]]
df.to_csv(sys.argv[3], sep = '\t', header = False, index = False)
