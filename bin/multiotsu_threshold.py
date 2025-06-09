# args: 'insert-sizes_chr22.txt'

import pandas as pd
import numpy as np
import sys
from skimage.filters import threshold_multiotsu

# file location for insert sizes
filepath = sys.argv[1]
# output thresholds location
outpath = sys.argv[2]
# read the insert sizes
df_22 = pd.read_csv(filepath, names = ['Length'])
# split into 4 classes: nfr, 1n, 2n, 3n
thresholds_22 = threshold_multiotsu(df_22['Length'].values, classes=4)
# save as df
pd.DataFrame(thresholds_22).to_csv(outpath, sep = '\t', header=False, index=False)