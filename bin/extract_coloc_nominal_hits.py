#!/usr/bin/env python3

import argparse
import os
import pandas as pd


NOMINAL_COLS = [
    "phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd", "n_var_in_cis", "dist_phe_var", "var_id",
    "var_chr", "var_from", "var_to", "nom_pval", "r_squared", "slope", "slope_se", "best_hit"
]


def norm_chr(val: str) -> str:
    s = str(val).strip()
    if s.lower().startswith("chr"):
        return s[3:]
    return s


def scan_nominal(path, keys_needed):
    out = {}
    if not os.path.exists(path):
        return out

    with open(path, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split()
            if len(parts) < 16:
                continue
            key = (parts[7], parts[0])  # (var_id, phe_id)
            if key in keys_needed and key not in out:
                out[key] = {
                    "slope": parts[13],
                    "slope_se": parts[14],
                    "nom_pval": parts[11],
                }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coloc", required=True)
    ap.add_argument("--nominal-dir-a", required=True)
    ap.add_argument("--nominal-dir-b", required=True)
    ap.add_argument("--out-nfr", required=True)
    ap.add_argument("--out-nuc", required=True)
    args = ap.parse_args()

    coloc = pd.read_csv(args.coloc, sep="\t", dtype=str)
    required = ["Chr", "hit1", "hit2", "NFR_peak", "Nucleosomal_peak"]
    missing = [c for c in required if c not in coloc.columns]
    if missing:
        raise ValueError(f"Missing required columns in coloc file: {missing}")

    coloc = coloc.copy()
    coloc["chr_num"] = coloc["Chr"].map(norm_chr)

    nfr_hits = {}
    nuc_hits = {}

    for chr_num, grp in coloc.groupby("chr_num", sort=False):
        nfr_keys = set(zip(grp["hit1"], grp["NFR_peak"]))
        nuc_keys = set(zip(grp["hit2"], grp["Nucleosomal_peak"]))

        nfr_file = os.path.join(args.nominal_dir_a, f"nominal_chr{chr_num}.txt")
        nuc_file = os.path.join(args.nominal_dir_b, f"nominal_chr{chr_num}.txt")

        for k, v in scan_nominal(nfr_file, nfr_keys).items():
            nfr_hits[(chr_num, k[0], k[1])] = v
        for k, v in scan_nominal(nuc_file, nuc_keys).items():
            nuc_hits[(chr_num, k[0], k[1])] = v

    nfr_rows = []
    nuc_rows = []

    for _, r in coloc.iterrows():
        chr_num = r["chr_num"]

        nfr_key = (chr_num, r["hit1"], r["NFR_peak"])
        nfr_match = nfr_hits.get(nfr_key, {"slope": None, "slope_se": None, "nom_pval": None})
        nfr_rows.append({
            "Chr": r["Chr"],
            "hit1": r["hit1"],
            "NFR_peak": r["NFR_peak"],
            "slope": nfr_match["slope"],
            "slope_se": nfr_match["slope_se"],
            "nom_pval": nfr_match["nom_pval"],
        })

        nuc_key = (chr_num, r["hit2"], r["Nucleosomal_peak"])
        nuc_match = nuc_hits.get(nuc_key, {"slope": None, "slope_se": None, "nom_pval": None})
        nuc_rows.append({
            "Chr": r["Chr"],
            "hit2": r["hit2"],
            "Nucleosomal_peak": r["Nucleosomal_peak"],
            "slope": nuc_match["slope"],
            "slope_se": nuc_match["slope_se"],
            "nom_pval": nuc_match["nom_pval"],
        })

    os.makedirs(os.path.dirname(args.out_nfr), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_nuc), exist_ok=True)

    pd.DataFrame(nfr_rows).to_csv(args.out_nfr, sep="\t", index=False)
    pd.DataFrame(nuc_rows).to_csv(args.out_nuc, sep="\t", index=False)


if __name__ == "__main__":
    main()
