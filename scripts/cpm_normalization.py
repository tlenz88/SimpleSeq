#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import os

df = pd.read_csv(sys.argv[1], sep = "\t", header = None)
mmr = df[2].sum() / 1e6
df[2] = df[2].div(mmr)
out = "".join([os.path.splitext(os.path.basename(sys.argv[1]))[0], "_cpm", os.path.splitext(os.path.basename(sys.argv[1]))[1]])
df.to_csv(out, sep = "\t", header = False, index = False)
