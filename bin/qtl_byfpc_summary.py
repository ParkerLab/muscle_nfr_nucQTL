#!/usr/bin/env python
# args: 1. celltype 2. nocov_dir 3. withcov_dir 4. outdir

# pkgs
import pandas as pd
import numpy as np
import csv
import os
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# inputs
celltype = sys.argv[1]

nocov_dir_nfr = sys.argv[2]
wcov_dir_nfr = sys.argv[3]
nocov_dir_nuc = sys.argv[4]
wcov_dir_nuc = sys.argv[5]
outdir = sys.argv[6]

fpc = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
def tidy_qtl(celltype, region, outdir, nocov_dir, wcov_dir):


fpc = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

    # track sigQTL numbers
    nocov_sigqtl = []
    wcov_sigqtl = []
    for i in fpc:
        fpc_num = i
        # get path
        cellcombo = celltype + '_' + region + '_F' + str(fpc_num) + '_G5/merged_QTLresults_fdrcorr.csv'
        nocov_results = nocov_dir + cellcombo
        wcov_results = wcov_dir + cellcombo
        # nocov results
        df1 = pd.read_csv(nocov_results, index_col=0).dropna()
        # find sig QTLs
        sig_df1 = df1[df1['adj_beta_qval'] < 0.05]
        nocov_sigqtl.append(len(sig_df1))
        # withcov results
        df2 = pd.read_csv(wcov_results, index_col=0).dropna()
        # find sig QTLs
        sig_df2 = df2[df2['adj_beta_qval'] < 0.05]
        wcov_sigqtl.append(len(sig_df2))
    # plot
    plt.plot(fpc, nocov_sigqtl, 'o-', label='Top5genoPCs + phenoPC + SNP')
    plt.plot(fpc, wcov_sigqtl, 'o-', label='Top5genoPCs + phenoPC + SNP + age + sex + BMI + batch + \n n_nuclei + TSS_enrichment + MedianFragLength_chr22')
    plt.title(celltype + ' ' + region)
    plt.ylabel('Number of Significant QTLs (FDR < 0.05)')
    plt.xlabel('Number of Phenotype/ATAC PCs')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 0.5))
    plt.savefig(outdir+celltype+"_"+region+"_summary.png",bbox_inches='tight')
    plt.close()
    return [nocov_sigqtl, wcov_sigqtl]


result_nfr = tidy_qtl(celltype, 'nfr', outdir, nocov_dir_nfr, wcov_dir_nfr)
result_nuc = tidy_qtl(celltype, 'nuc', outdir, nocov_dir_nuc, wcov_dir_nuc)


# save a summary df
summary_df = pd.DataFrame({celltype+'_num_sigQTL_nfr_nocov': result_nfr[0],
                           celltype+'_num_sigQTL_nfr_wcov': result_nfr[1],
                           celltype+'_num_sigQTL_nuc_nocov': result_nuc[0],
                           celltype+'_num_sigQTL_nuc_wcov': result_nuc[1]
                          },
                         index = fpc)
summary_df.to_csv(outdir+celltype+'_qtl_by_fpcnum_summary.txt', sep = '\t', index_label='FPC')

