#!/usr/bin/env python3
"""
Filter a featureCounts-style peak x sample matrix for eligible donors.
Modes:
1) split mode: random A/B split (legacy)
2) single-output mode: write one filtered matrix
   - supports deterministic batch-based subset selection
"""

import argparse
import sys
import pandas as pd
import numpy as np


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--raw-counts", required=True)
    p.add_argument("--sample-info", required=True)
    p.add_argument("--keep-peaks", default=None)
    p.add_argument("--keep-peaks-col", default=None)
    p.add_argument("--out", default=None)
    p.add_argument("--out-splitA", default=None)
    p.add_argument("--out-splitB", default=None)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--exclude-donors", nargs="*", default=[])
    p.add_argument("--subsetA-batches", nargs="*", default=None)
    p.add_argument("--subset-label", choices=["A", "B"], default=None)
    return p.parse_args()


def load_keep_peaks(path, col=None):
    if path is None:
        return None
    if col:
        kp_df = pd.read_csv(path, sep="\t", dtype=str)
        if col not in kp_df.columns:
            raise ValueError(f"keep_peaks_col not found: {col}")
        vals = kp_df[col].dropna().astype(str).tolist()
        vals = [x for x in vals if x and x != "nan"]
    else:
        with open(path, "r") as fh:
            vals = [line.strip() for line in fh if line.strip()]
    return set(vals)


def main():
    args = parse_args()

    do_split = args.out_splitA is not None and args.out_splitB is not None
    do_single = args.out is not None

    if do_split == do_single:
        raise ValueError("Provide either --out or --out-splitA/--out-splitB")

    sample_info = pd.read_csv(args.sample_info, sep="\t", dtype=str)
    if "SAMPLE" not in sample_info.columns:
        raise ValueError("sample_info missing required column: SAMPLE")
    if "batch" not in sample_info.columns:
        raise ValueError("sample_info missing required column: batch")

    sample_info = sample_info.copy()
    sample_info["SAMPLE"] = sample_info["SAMPLE"].astype(str)
    sample_info["batch"] = sample_info["batch"].astype(str)

    df = pd.read_csv(args.raw_counts, sep="\t", dtype=str)
    if df.shape[1] < 2:
        raise ValueError("raw counts must have at least 2 columns")

    peak_col = df.columns[0]
    sample_cols = list(df.columns[1:])

    # normalization.R expects donor IDs in headers; strip suffix after "__"
    col_to_donor = {c: c.split("__", 1)[0] for c in sample_cols}

    exclude = set(str(x) for x in args.exclude_donors)
    donors_in_info = set(sample_info["SAMPLE"].tolist())
    donors_in_counts = set(col_to_donor.values())
    eligible_donors = sorted(
        d for d in donors_in_counts if d in donors_in_info and d not in exclude
    )

    if not eligible_donors:
        raise ValueError("No eligible donors after intersection/exclusion")

    keep_set = load_keep_peaks(args.keep_peaks, args.keep_peaks_col)
    if keep_set is not None:
        df = df[df[peak_col].isin(keep_set)]
        present_peaks = set(df[peak_col].tolist())
        missing = sorted(keep_set - present_peaks)
        if missing:
            print(
                f"WARNING: {len(missing)} peaks from keep_peaks missing in counts",
                file=sys.stderr,
            )
            print("Missing (first 10): " + ", ".join(missing[:10]), file=sys.stderr)

    donor_to_batch = sample_info.set_index("SAMPLE")["batch"].to_dict()

    if do_split:
        rng = np.random.RandomState(args.seed)
        donors_A = set()
        donors_B = set()

        batch_to_donors = {}
        for d in eligible_donors:
            b = donor_to_batch.get(d, "NA")
            batch_to_donors.setdefault(b, []).append(d)

        for donors in batch_to_donors.values():
            donors = sorted(donors)
            rng.shuffle(donors)
            nA = len(donors) // 2
            donors_A.update(donors[:nA])
            donors_B.update(donors[nA:])

        cols_A = [c for c in sample_cols if col_to_donor[c] in donors_A]
        cols_B = [c for c in sample_cols if col_to_donor[c] in donors_B]

        if not cols_A or not cols_B:
            raise ValueError("Split resulted in empty A or B column set")

        outA = df[[peak_col] + cols_A].copy()
        outA.columns = [peak_col] + [col_to_donor[c] for c in cols_A]
        outA.to_csv(args.out_splitA, sep="\t", index=False)

        outB = df[[peak_col] + cols_B].copy()
        outB.columns = [peak_col] + [col_to_donor[c] for c in cols_B]
        outB.to_csv(args.out_splitB, sep="\t", index=False)

        print(f"Eligible donors: {len(eligible_donors)}", file=sys.stderr)
        print(f"SplitA columns: {len(cols_A)}", file=sys.stderr)
        print(f"SplitB columns: {len(cols_B)}", file=sys.stderr)
        return

    # single-output mode
    if args.subsetA_batches is not None and args.subset_label is None:
        raise ValueError("--subset-label is required when --subsetA-batches is used")
    if args.subset_label is not None and args.subsetA_batches is None:
        raise ValueError("--subsetA-batches is required when --subset-label is used")

    if args.subsetA_batches is not None:
        subsetA_batches = set(str(x) for x in args.subsetA_batches)
        if args.subset_label == "A":
            selected_donors = {
                d for d in eligible_donors if donor_to_batch.get(d, "") in subsetA_batches
            }
        else:
            selected_donors = {
                d for d in eligible_donors if donor_to_batch.get(d, "") not in subsetA_batches
            }
    else:
        selected_donors = set(eligible_donors)

    cols_keep = [c for c in sample_cols if col_to_donor[c] in selected_donors]
    if not cols_keep:
        raise ValueError("No eligible columns after filtering")

    out = df[[peak_col] + cols_keep].copy()
    out.columns = [peak_col] + [col_to_donor[c] for c in cols_keep]
    out.to_csv(args.out, sep="\t", index=False)

    print(f"Eligible donors: {len(eligible_donors)}", file=sys.stderr)
    print(f"Selected donors: {len(selected_donors)}", file=sys.stderr)
    print(f"Columns kept: {len(cols_keep)}", file=sys.stderr)


if __name__ == "__main__":
    main()
