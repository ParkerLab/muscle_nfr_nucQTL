# args: 1. input_count_mat 2. threshold_score 3. threshold_percentage 4.output_dir

# import pkgs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyreadr
import sys

# inputs
input_count_mat = sys.argv[1]
threshold_score = float(sys.argv[2])
threshold_percentage = float(sys.argv[3])
output = sys.argv[4]

# read in
df_nucleo_x = pyreadr.read_r(input_count_mat)[None]

def filter_peaks(df, threshold_score, threshold_percentage):
    # convert raw data to avg df
    #avg_df = df.set_index('name')
    # length as 151bp
    #avg_df = avg_df.iloc[:,2:]/101
    # init for qualified peaks
    qualified_peaks = []
    # threshold
    #threshold_score = 1.0
    threshold_count = int(df.shape[1]*threshold_percentage)
    for i in range(df.shape[0]):
        # check nuc peak one by one
        peak_occ_list = list(df.iloc[i])
        peak_name = df.index[i]
        if sum(i >= threshold_score for i in peak_occ_list) >= threshold_count:
            qualified_peaks.append(peak_name)
    return qualified_peaks

qualified_peaks = filter_peaks(df_nucleo_x, threshold_score, threshold_percentage)

df_nucleo_x.loc[qualified_peaks].to_csv(output, sep = '\t')