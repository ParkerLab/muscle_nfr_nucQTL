# pkgs
import pandas as pd
import csv
import os
import sys

# inputs
summary = sys.argv[1]
celltype = sys.argv[2]
#nocov_results_dir = '../../results/results_05082024/QTLscan_nocov_permutation/'
output = sys.argv[3]

# read in the df
summary_df = pd.read_csv(summary, sep = '\t', header=0)
# get the col name - looking at nocov results
nfr_col = celltype + '_num_sigQTL_nfr_nocov'
nuc_col = celltype + '_num_sigQTL_nuc_nocov'
# get max fpc
nfr_fpc = summary_df['FPC'][summary_df[nfr_col].idxmax()]
nuc_fpc = summary_df['FPC'][summary_df[nuc_col].idxmax()]
# output

d = {'nfr': [nfr_fpc], 'nuc': [nuc_fpc]}
df = pd.DataFrame(data=d, index=[celltype])

df.to_csv(output)