import pandas as pd
import csv
import numpy as np
import glob
import sys

# parameters
celltype = sys.argv[1]
sample_input = sys.argv[2]
outdir = sys.argv[3]

pattern = "../../data/bams-by-cluster-sample/" + celltype
pattern += "/chr22/sizes/atac-"
pattern += celltype
pattern += "@fusion@*.txt"
filelist = glob.glob(pattern)

# read sample info
sample_info = pd.read_csv(sample_input, sep = '\t')
celltype_sample_info = sample_info[sample_info['coarse_cluster_name'] == 'fusion.'+celltype]

# size df
size_df = pd.DataFrame(data=filelist, columns = ['path'])
size_df['SNG.1ST'] = [i.split('@fusion@')[1].split('--')[0] for i in size_df['path'].tolist()]

medians = []
for i in size_df['path']:
    #print(i)
    chr22_sizes = pd.read_csv(i, names=['size'])['size']
    medians.append(np.median(chr22_sizes))

size_df['chr22_median_len'] = medians
celltype_sample_info = celltype_sample_info.merge(size_df)
celltype_sample_info.to_csv(outdir, sep = '\t')