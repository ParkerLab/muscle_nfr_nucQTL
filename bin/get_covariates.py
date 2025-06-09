# pkgs
import pandas as pd
import csv
import numpy as np
import pyreadr
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import sys

fpc_num, gpc_num, outdir, figdir = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
sample_info = sys.argv[5]
vcf_eigenvec = sys.argv[6]
counts = sys.argv[7]
celltype = sys.argv[8]
region = sys.argv[9]
# outputs files
outdir1 = outdir+"qtl_cov_full_"+celltype+"_"+region+"_FPC" + str(fpc_num) + "_GPC" + str(gpc_num) + ".txt"
outdir2 = outdir+"qtl_cov_simple_"+celltype+"_"+region+"_FPC" + str(fpc_num) + "_GPC" + str(gpc_num) + ".txt"
# read in
df_nfr_x = pyreadr.read_r(counts)[None].T
#df_nfr_x = pyreadr.read_r('../../occ_matrix_nucmap_50ext_1-092622/filtered_occ_matrix_nucmap_50ext_1_09252022_normalized.rds')[None]
# pca analysis
pca = PCA(n_components=fpc_num)
pca.fit(df_nfr_x)

# Reformat and view results
loadings = pd.DataFrame(pca.components_.T,
columns=['FPC%s' % i for i in range(fpc_num)], index=df_nfr_x.columns)
#loadings['batch'] = [i.split('_')[2] for i in loadings.index]
#loadings = loadings.rename(index=lambda s: s.split('__')[0])

# bring in sample df
sample_df_newinfo = pd.read_csv(sample_info, sep ='\t')
sample_df_newinfo = sample_df_newinfo[['SAMPLE', 'batch', 'age', 'sex', 'bmi', 'n_nuclei', 'tss_enrichment', 'chr22_median_len']]
sample_df_newinfo = sample_df_newinfo.drop_duplicates()

# create lkey for loadings
loadings['sample'] = loadings.index
# merge loadings and sample df
loadings = loadings.merge(sample_df_newinfo, left_on='sample', right_on='SAMPLE')
# set sample id as index
loadings = loadings.set_index('sample')
loadings.index.name = None

# adjust covariates
# convert sex to numeric
loadings['sex'].replace(['F','M'], [0,1], inplace=True)
# same output
# spread the batch info
loadings = pd.get_dummies(loadings, prefix=['batch'], columns=['batch'], drop_first=True, dtype='int')

# add genotype PCs
# reference from the plink_PCA_analysis.ipynb
vec = pd.read_csv(vcf_eigenvec, sep = '\t')
# specify genotype PC col names
vec.columns = vec.columns.str.replace('PC', 'GPC')
# select number of genoPCs based on input
vec = pd.concat([vec.iloc[:,0:2], vec.iloc[:,2:2+gpc_num]], axis=1)
vec = vec.merge(loadings, left_on='IID', right_on='SAMPLE', how = 'inner')
# set index
# set sample id as index
vec = vec.set_index('IID')
vec.index.name = 'id'
# remove SAMPLE
vec = vec.drop(['SAMPLE', '#FID'], axis=1)
vec = vec.astype(float)
# add notes to outdir
#outdir = outdir+'/'+inputdir.split('/')[-1].split('.')[0] + '_cov.txt'
# first version
vec.T.to_csv(outdir1,sep = '\t', index_label='id')
# drop other covs
nocov_df = vec
# drop covs
nocov_df = nocov_df.drop(columns = ['sex','age','bmi', 'n_nuclei', 'tss_enrichment', 'chr22_median_len'])
# drop batch
nocov_df = nocov_df.loc[:, ~nocov_df.columns.str.startswith('batch')]
# save
nocov_df.T.to_csv(outdir2,sep = '\t', index_label='id')

# read cov files
cov_df = vec.T
cov_corr = cov_df.T.iloc[:, gpc_num:gpc_num+fpc_num+15].corr(method = 'spearman')
# select cols
cov_corr = cov_corr[['age', 'sex', 'bmi', 'n_nuclei', 'tss_enrichment', 'chr22_median_len',
                     'batch_1', 'batch_2', 'batch_3', 'batch_4', 
                     'batch_5', 'batch_6', 'batch_7', 'batch_8', 'batch_9']]
# select rows
cov_corr = cov_corr.iloc[0:fpc_num]
# plot heatmap
plt.figure(figsize=(12,12))
plt.title(celltype+"_"+region)
sns.heatmap(cov_corr, cmap='vlag', annot=True, vmin=-1, vmax=1)
#plt.show()
fname = figdir + 'FPC' + str(fpc_num) + '_' + celltype + '_' + region + '.png'
plt.savefig(fname)