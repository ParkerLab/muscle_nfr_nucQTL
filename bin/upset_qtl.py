# args: 

# import pkgs
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import upsetplot
from upsetplot import from_memberships
from upsetplot import plot
from upsetplot import UpSet
from matplotlib import pyplot
from matplotlib_venn import venn2
import sys
from matplotlib import font_manager
# process input
ca_nfr_path = sys.argv[1]
ca_nuc_path = sys.argv[2]
nfr_nuc_path = sys.argv[3]
nfr_ca_path=sys.argv[4]
nuc_ca_path=sys.argv[5]
nuc_nfr_path=sys.argv[6]
all_path=sys.argv[7]
celltype = sys.argv[8]
nfr = sys.argv[9]
nuc = sys.argv[10]
ca = sys.argv[11]
output = sys.argv[12]
# read in data
nfr_ca_df = pd.read_csv(nfr_ca_path, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])
nuc_ca_df = pd.read_csv(nuc_ca_path, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])
nfr_nuc_df = pd.read_csv(nfr_nuc_path, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])

ca_nfr_df = pd.read_csv(ca_nfr_path, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])
ca_nuc_df = pd.read_csv(ca_nuc_path, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])
nuc_nfr_df = pd.read_csv(nuc_nfr_path, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])

all_df = pd.read_csv(all_path, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])

caqtl_df = pd.read_csv(ca, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])
nfrqtl_df = pd.read_csv(nfr, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])
nucqtl_df = pd.read_csv(nuc, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])
# Raise warning if length does not equal
if len(nfr_ca_df) != len(ca_nfr_df):
    print("nfr_ca does not equal ca_nfr")
if len(nuc_ca_df) != len(ca_nuc_df):
    print("nuc_ca does not equal ca_nuc")
if len(nfr_nuc_df) != len(nuc_nfr_df):
    print("nfr_nuc does not equal nuc_nfr")
# get the rest
ca_only = len(caqtl_df) - len(ca_nfr_df) - len(ca_nuc_df) + len(all_df)
nfr_only = len(nfrqtl_df) - len(nfr_ca_df) - len(nfr_nuc_df) + len(all_df)
nuc_only = len(nucqtl_df) - len(nuc_nfr_df) - len(nuc_ca_df) + len(all_df)
data = {
    'caQTL': [1, 0, 0, 1, 1, 0, 1],
    'nfrQTL': [0, 1, 0, 1, 0, 1, 1],
    'nucQTL': [0, 0, 1, 0, 1, 1, 1],
}
# add counts
counts = [ca_only, nfr_only, nuc_only,
          len(nfr_ca_df)-len(all_df), len(ca_nuc_df)-len(all_df), len(nfr_nuc_df)-len(all_df),
          len(all_df)]
# Create DataFrame
df_combinations = pd.DataFrame(data)

# Create a Series with a MultiIndex based on the combinations
index = pd.MultiIndex.from_frame(df_combinations)
series = pd.Series(counts, index=index)
# plot
font_path = '/home/xiaoouw/font/Arial.ttf'
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)

plt.rcParams['font.family'] = 'Arial'

fig = plt.gcf()  # Get the current figure
fig.set_size_inches(2, 3)

ax_dict = UpSet(series, subset_size="sum", show_counts='{:,.0f}', sort_by='cardinality').plot()
# Add a title to the plot
plt.suptitle(celltype, fontsize=12)
pyplot.savefig(output, dpi=200)