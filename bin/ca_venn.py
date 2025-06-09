# args: 1. ca nfr bed path 2. ca nuc bed path 3. figsave location 4. celltype 5. caQTL

# import pkgs
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import sys

# process input
ca_nfr_path = sys.argv[1]
ca_nuc_path = sys.argv[2]
output = sys.argv[3]
celltype = sys.argv[4]
caqtl = sys.argv[5]

caqtl_df = pd.read_csv(caqtl, sep = '\t')
num_caqtl = len(caqtl_df )

ca_nfr_df = pd.read_csv(ca_nfr_path, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])
ca_nuc_df = pd.read_csv(ca_nuc_path, sep = '\t', names=['phe_chr', 'phe_from', 'phe_to'])

# Convert dataframe rows to sets of tuples
nfr_tuples = set(ca_nfr_df.apply(tuple, axis=1))
nuc_tuples = set(ca_nuc_df.apply(tuple, axis=1))

# Create the venn diagram
plt.figure(figsize=(6, 6))
venn = venn2([nfr_tuples, nuc_tuples], ('nfr', 'nuc'))

plt.title("Venn Diagram of " + celltype + " ca partition, #caQTL=" + str(num_caqtl))
plt.show()
# save
plt.savefig(output)