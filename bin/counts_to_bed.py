# import pkgs
import pandas as pd
import csv
import numpy as np
import pyreadr
import sys

# load inputs and outputs
inputdir = sys.argv[1]
outdir = sys.argv[2]

# read in
df_nfr_x = pyreadr.read_r(inputdir)[None].T
# rename columns
df_nfr_x = df_nfr_x.rename(columns=lambda s: s.split('__')[0])
# add the 6 cols for bed file
df_nfr_x['#Chr'] = [i.split('.')[0] for i in list(df_nfr_x.index)]
df_nfr_x['start'] = [i.split('.')[1] for i in list(df_nfr_x.index)]
# mistake here! should be [2]!!!
#df_nfr_x['end'] = [i.split('.')[1] for i in list(df_nfr_x.index)]
df_nfr_x['end'] = [i.split('.')[2] for i in list(df_nfr_x.index)]
df_nfr_x['pid'] = [i for i in list(df_nfr_x.index)]
df_nfr_x['gid'] = [i for i in list(df_nfr_x.index)]
df_nfr_x['strand'] = '.'
# rearrange columns
cols = df_nfr_x.columns.tolist()
cols = cols[-6:] + cols[:-6]
df_nfr_x = df_nfr_x[cols]
# sort the df
df_nfr_x['start'] = df_nfr_x['start'].astype(int)
df_nfr_x['end'] = df_nfr_x['end'].astype(int)
df_nfr_x = df_nfr_x.sort_values(by=['#Chr', 'start'])
# split by chr
for i in range(22):
    # chr number
    j = i+1
    chr_string = 'chr' + str(j)
    df = df_nfr_x[df_nfr_x['#Chr'] == chr_string]
    # save chr bed
    #df.to_csv(outdir+chr_string+'.bed', sep = '\t', header = True, index = None)
# save as bed file
df_nfr_x.to_csv(outdir, sep = '\t', header = True, index = None)