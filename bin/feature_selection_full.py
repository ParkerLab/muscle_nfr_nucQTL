# args: 1. input_count_mat 2. nfrPeaks file 3.output_dir

# import pkgs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyreadr
import sys

# inputs
input_count_mat = sys.argv[1]
nfrPeaks_dir = sys.argv[2]
output = sys.argv[3]

# read in
df_nucleo_x = pyreadr.read_r(input_count_mat)[None]

nfr_df = pd.read_csv(nfrPeaks_dir, sep = '\t', index_col=0)

qualified_peaks = nfr_df.index.tolist()

df_nucleo_x.loc[qualified_peaks].to_csv(output, sep = '\t')