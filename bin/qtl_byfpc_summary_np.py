#!/usr/bin/env python
"""
Summarize the number of significant caQTLs across FPC settings for narrowPeak runs.

Args:
  1) celltype
  2) qtl_dir (directory containing {celltype}_ca_F{fpc}_G{gpc}/merged_QTLresults_fdrcorr.csv)
  3) outdir
  4) gpc (default "5" if omitted)
"""

import os
import sys
import pandas as pd


FPC_VALUES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]


def count_sig_qtl(result_file: str, qval_col: str = "adj_beta_qval", threshold: float = 0.05) -> int:
    df = pd.read_csv(result_file, index_col=0).dropna()
    if qval_col not in df.columns:
        raise ValueError(f"Column '{qval_col}' not found in {result_file}")
    return int((df[qval_col] < threshold).sum())


def main() -> None:
    if len(sys.argv) < 4:
        raise SystemExit("Usage: qtl_byfpc_summary_np.py <celltype> <qtl_dir> <outdir> [gpc]")

    celltype = sys.argv[1]
    qtl_dir = sys.argv[2]
    outdir = sys.argv[3]
    gpc = sys.argv[4] if len(sys.argv) > 4 else "5"

    sig_counts = []
    for fpc in FPC_VALUES:
        result_file = os.path.join(
            qtl_dir,
            f"{celltype}_ca_F{fpc}_G{gpc}",
            "merged_QTLresults_fdrcorr.csv",
        )
        sig_counts.append(count_sig_qtl(result_file))

    # Keep legacy column name for compatibility with find_max_caqtl.py.
    summary_df = pd.DataFrame(
        {
            f"{celltype}_num_sigQTL_nfr_nocov_75ext": sig_counts,
            f"{celltype}_num_sigQTL_ca_nocov_narrowpeak": sig_counts,
        },
        index=FPC_VALUES,
    )

    os.makedirs(outdir, exist_ok=True)
    summary_df.to_csv(
        os.path.join(outdir, f"{celltype}_qtl_by_fpcnum_summary.txt"),
        sep="\t",
        index_label="FPC",
    )


if __name__ == "__main__":
    main()
