#!/usr/bin/env python3

import sys
import pandas as pd
import pyreadr
from sklearn.decomposition import PCA


# args:
# 1 fpc_num
# 2 gpc_num
# 3 output_simple_cov
# 4 vcf_eigenvec
# 5 normalized_counts_rds
fpc_num = int(sys.argv[1])
gpc_num = int(sys.argv[2])
out_simple = sys.argv[3]
vcf_eigenvec = sys.argv[4]
counts_rds = sys.argv[5]

# normalized matrix from normalization.R is sample x peak; transpose for PCA-on-features pattern
mat = pyreadr.read_r(counts_rds)[None].T

pca = PCA(n_components=fpc_num)
pca.fit(mat)

loadings = pd.DataFrame(
    pca.components_.T,
    columns=[f"FPC{i}" for i in range(fpc_num)],
    index=mat.columns,
)
loadings["sample"] = loadings.index

vec = pd.read_csv(vcf_eigenvec, sep="\t")
vec.columns = vec.columns.str.replace("PC", "GPC")
vec = pd.concat([vec.iloc[:, 0:2], vec.iloc[:, 2 : 2 + gpc_num]], axis=1)

merged = vec.merge(loadings, left_on="IID", right_on="sample", how="inner")
merged = merged.set_index("IID")
merged.index.name = "id"

drop_cols = [c for c in ["sample", "#FID"] if c in merged.columns]
if drop_cols:
    merged = merged.drop(columns=drop_cols)

# simple cov = genotype PCs + feature PCs only
keep_cols = [c for c in merged.columns if c.startswith("GPC") or c.startswith("FPC")]
simple = merged[keep_cols].astype(float)

simple.T.to_csv(out_simple, sep="\t", index_label="id")
