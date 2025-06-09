#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas
import numpy
import math
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
import glob
import time
import pybedtools
import subprocess as sp
dpi = 150
matplotlib.rcParams['figure.dpi']= dpi
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams.update({'font.size': 12})
#import pandas_extra
#from pandas_extras import ExtraFunctions
#from plot_utils import make_raincloud
from scipy import stats
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu
import random
import re
import argparse


def getOpts():
    parser = argparse.ArgumentParser(description='plot LDSC joint enrichment')
    parser.add_argument('--joint', required=True, nargs='+', help ="""ldsc outputs """)    
    args = parser.parse_args()
    return args


# In[2]:
args = getOpts()



# In[2]:


"""Fix the directory where all figures will be saved"""

fdir = "."    
#pe = pandas_extra.ExtraFunctions(fdir)
pe = fdir

# ukbb info
#u = pandas.read_csv("/lab/work/arushiv/muscle-sn/data/ldsc/ukbb/data/LDSC Sumstat Manifest for Neale UKB GWAS - ukb31063_ldsc_sumstat_manifest.tsv",
#                    sep='\t')
u = pandas.read_csv("/home/xiaoouw/Muscle_snATAC_nucQTL/scripts/nextflow/ldsc_allgwas_consensus-cluster/data/LDSC Sumstat Manifest for Neale UKB GWAS - ukb31063_ldsc_sumstat_manifest.tsv",
                    sep='\t')
#print(u[u['ldsc_sumstat_file'].str.contains("20502")])
                    #usecols=['phenotype', 'description'], sep='\t').drop_duplicates()
u = u[(u['is_primary_gwas']) & (u['ldsc_confidence']=="high") & (u['ldsc_h2_significance'] == "z7")][['ldsc_sumstat_file', 'description']]
u = u[u['ldsc_sumstat_file'].str.contains("both_sexes")]
#u['trait'] = u['ldsc_sumstat_file'].str.replace(".bgz", "").str.replace(".gz", "").str.replace(".ldsc.imputed_v3.both_sexes.tsv", "").replace(".imputed_v3.ldsc.both_sexes.tsv", "")
u['trait'] = u['ldsc_sumstat_file'].map(lambda x: os.path.basename(x).split(".")[0])
u = u[['trait', 'description']].drop_duplicates()
# shorten some labels:
u['description'] = u['description'].str.replace(", not elsewhere classified", "")
u['description'] = u['description'].map(lambda x: re.sub(r'.*doctor: ', '', x))
u = u[~ u['description'].str.contains("None of the above")]
u = u[~u['description'].duplicated()]
print(u.shape)
#print(u[u['trait'].str.contains("20502")])
u.head()


# In[13]:


u.to_csv("gwas_shortlist_for_indexing.tsv", sep='\t', index=False)


# In[14]:

# joint enrich
#annot_string = "ld_scores_by_annot/Adipocyte.annot.@.ld,ld_scores_by_annot/Adipose-bulk.annot.@.ld,ld_scores_by_annot/Endothelial.annot.@.ld,ld_scores_by_annot/Islet-beta.annot.@.ld,ld_scores_by_annot/Macrophage.annot.@.ld,ld_scores_by_annot/Mesenchymal_Stem_Cell.annot.@.ld,ld_scores_by_annot/Muscle_fiber.annot.@.ld,ld_scores_by_annot/Neuronal.annot.@.ld,ld_scores_by_annot/Satellite_Cell.annot.@.ld,ld_scores_by_annot/Smooth_Muscle.annot.@.ld,ld_scores_by_annot/T_cell.annot.@.ld"
annot_string = "ld_scores_by_annot/NucleoATAC_NFR_sigQTL.annot.@.ld, ld_scores_by_annot/NucleoATAC_Nucleo_sigQTL.annot.@.ld"
annots = annot_string.replace("ld_scores_by_annot/", "").replace(".annot.@.ld", "").split(',')


# In[15]:


dlist = []
for fname in args.joint:
    
    f = pandas.read_csv(fname, sep='\t')
    trait = os.path.basename(fname).replace(".joint.results", "")
    f['trait'] = trait
    f = f[f['Category'].str.startswith("L2_")]

    d = {}
    for i, x in enumerate(annots):
        d[f"L2_{i}"] = x

    f['cluster'] = f['Category'].map(d)

    #f.to_csv(f"{trait}.joint.tsv", sep='\t', index=False)
    dlist.append(f)
    
d = pandas.concat(dlist)
d = pandas.merge(d, u, how="inner", on="trait")
d['fdr_by_p'] = multipletests(d['Enrichment_p'], method="fdr_by")[1]
d.to_csv("joint_enrichment.tsv", sep='\t', na_rep="NA")

d.head()

# Select traits that have at least one annot significant in the joint model, and enrichment > 10
selected_traits = d[((d['Enrichment']>10) | (d['Enrichment']<0)) & (d['fdr_by_p']<0.05) & (((d['Coefficient']>0) & ((d['Coefficient'] - 1.96*d['Coefficient_std_error']) > 0)) | ((d['Coefficient']<0) & ((d['Coefficient'] + 1.96*d['Coefficient_std_error']) < 0)))]['description'].drop_duplicates().tolist()

print(f"selecting {len(selected_traits)} traits")
d = d[d['description'].isin(selected_traits)]

selected_clusters = annots

# In[16]:


d1 = d.pivot_table(index="cluster", columns="description", values="Coefficient_z-score")
d1.to_csv("ldsc_joint_coeff_fdr_by_p.tsv", sep='\t', na_rep="NA")

dna = d[(((d['Coefficient']>0) & ((d['Coefficient'] - 1.96*d['Coefficient_std_error']) > 0)) | ((d['Coefficient']<0) & ((d['Coefficient'] + 1.96*d['Coefficient_std_error']) < 0)))].pivot_table(index="cluster", columns="description", values="Coefficient_z-score").fillna(0)
dna.to_csv("ldsc_joint_coeff_fdr_by_p_nas.tsv", sep='\t', na_rep="NA")

d1.head()


# In[7]:


g = sns.clustermap(d1, cmap="viridis", figsize=(40,12))
pe.save("fig.significant_coeff_joint_1.png")

row_order = [i.get_text() for i in g.ax_heatmap.yaxis.get_ticklabels()]
col_order = [i.get_text() for i in g.ax_heatmap.xaxis.get_ticklabels()]


# In[8]:


maskvals = d[((d['Coefficient']>0) & ((d['Coefficient'] - 1.96*d['Coefficient_std_error']) > 0)) | ((d['Coefficient']<0) & ((d['Coefficient'] + 1.96*d['Coefficient_std_error']) < 0))]
d3 = maskvals.pivot_table(index="cluster", columns="description", values="Coefficient").fillna("")
d3 = d3.applymap(lambda x: "*" if x!="" else x)
selected_traits = d3.columns
selected_clusters = d3.index
d3


# In[9]:


d1 = d[(d['cluster'].isin(selected_clusters)) & (d['description'].isin(selected_traits))]
d2 = d1.pivot_table(index="cluster", columns="description", values="Coefficient_z-score")
d2


# In[10]:


print(d2.shape)
print(d3.shape)


# In[11]:


plt.figure(figsize=(30, 5))
sns.heatmap(d2, cmap="vlag", center=0, annot=d3, fmt="", annot_kws={'size': 10})
pe.saveb("fig.gwas_coeff_z_scores_sig-conf-intervals.png")


# In[ ]:




