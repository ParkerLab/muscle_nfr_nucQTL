# pkgs
import pandas as pd
import numpy as np
import csv
import os
#import statsmodels.stats.multitest as smm
import matplotlib.pyplot as plt
import sys

# read the file
directory=sys.argv[1]
pvalue=float(sys.argv[2])
# init empty df
headers = 'phe_id phe_chr phe_from phe_to phe_strd n_var_in_cis dist_phe_var var_id var_chr var_from var_to dof1 dof2 bml1 bml2 nom_pval r_squared slope slope_se adj_emp_pval adj_beta_pval'
df = pd.DataFrame(columns=headers.split(" "))
# read the threshold files
#directory = '../QTLresults/NFR_multiotsu' # change to relative path if necessary
# iterate over files in
# that directory
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    # read the chr df
    if os.path.isfile(f) and 'chr' in f:
        # note the df is space-separated
        chr_df = pd.read_csv(f, sep = ' ', names = headers.split(" "))
        # append the chr df
        df = pd.concat([df, chr_df])
# make sure certain columns are number
df['adj_emp_pval'] = df['adj_emp_pval'].astype(float)
df['adj_beta_pval'] = df['adj_beta_pval'].astype(float)
df['phe_from'] = df['phe_from'].astype(int)
df['phe_to'] = df['phe_to'].astype(int)
# sort the df
df = df.sort_values(by=['phe_chr', 'phe_from','phe_to'])
# save the file
df.to_csv(directory+'/merged_QTLresults.txt', sep ='\t', header = True, index = None)
# report some outputs
#print('Total number of peaks: '+str(len(df)))
#print('Filtering by p-values from permutation: '+str(len(df[df['adj_emp_pval'] <= pvalue])))
#print('Filtering by p-values given by the fitted beta distribution: '+str(len(df[df['adj_beta_pval'] < pvalue])))
#print('Filtering by p-values from both columns: '+str(len(df[(df['adj_emp_pval'] < pvalue) & (df['adj_beta_pval'] < pvalue)])))