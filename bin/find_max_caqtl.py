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
ca_col = celltype + '_num_sigQTL_nfr_nocov_75ext'
# get max fpc
ca_fpc = summary_df['FPC'][summary_df[ca_col].idxmax()]
# output

d = {'ca': [ca_fpc]}
df = pd.DataFrame(data=d, index=[celltype])

df.to_csv(output)