#!/usr/bin/env python3

import argparse
import gzip
import os
import pandas as pd


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


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info_field(info, key):
    prefix = key + "="
    for item in info.split(";"):
        if item.startswith(prefix):
            return item[len(prefix):]
    return None


def first_float(x):
    if x is None:
        return None
    if "," in x:
        x = x.split(",")[0]
    try:
        return float(x)
    except Exception:
        return None


def annotate_from_vcf(vcf_path, targets):
    ann = {}

    for v in targets:
        parts = v.split(":")
        if len(parts) == 4:
            ann[v] = {
                "variant_id": v,
                "var_chr": parts[0],
                "var_pos": parts[1],
                "other_allele": parts[2],
                "effect_allele": parts[3],
                "eaf": None,
            }

    need = {v for v in targets if v not in ann}

    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, _qual, _filt, info = cols[:8]
            alts = alt.split(",")
            af_val = parse_info_field(info, "AF")
            af_first = first_float(af_val)

            if vid and vid != ".":
                for one_id in vid.split(";"):
                    if one_id in need and one_id not in ann:
                        ann[one_id] = {
                            "variant_id": one_id,
                            "var_chr": chrom,
                            "var_pos": pos,
                            "other_allele": ref,
                            "effect_allele": alts[0],
                            "eaf": af_first,
                        }

            for i, a in enumerate(alts):
                key = f"{chrom}:{pos}:{ref}:{a}"
                if key in targets:
                    af_i = None
                    if af_val is not None:
                        parts = af_val.split(",")
                        if i < len(parts):
                            af_i = first_float(parts[i])
                    ann[key] = {
                        "variant_id": key,
                        "var_chr": chrom,
                        "var_pos": pos,
                        "other_allele": ref,
                        "effect_allele": a,
                        "eaf": af_i,
                    }

    return ann


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coloc", required=True)
    ap.add_argument("--nominal-dir-a", required=True, help="NFR nominal directory (splitA)")
    ap.add_argument("--nominal-dir-b", required=True, help="NUC nominal directory (splitB)")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--n-exp", required=True, type=int)
    ap.add_argument("--n-out", required=True, type=int)
    ap.add_argument("--out", required=True)
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
        nfr_keys = set(zip(grp["hit1"], grp["NFR_peak"])) | set(zip(grp["hit2"], grp["NFR_peak"]))
        nuc_keys = set(zip(grp["hit1"], grp["Nucleosomal_peak"])) | set(zip(grp["hit2"], grp["Nucleosomal_peak"]))

        nfr_file = os.path.join(args.nominal_dir_a, f"nominal_chr{chr_num}.txt")
        nuc_file = os.path.join(args.nominal_dir_b, f"nominal_chr{chr_num}.txt")

        nfr_hits.update({(chr_num, k[0], k[1]): v for k, v in scan_nominal(nfr_file, nfr_keys).items()})
        nuc_hits.update({(chr_num, k[0], k[1]): v for k, v in scan_nominal(nuc_file, nuc_keys).items()})

    targets = set(coloc["hit1"].dropna().astype(str)) | set(coloc["hit2"].dropna().astype(str))
    ann = annotate_from_vcf(args.vcf, targets)

    rows = []
    for _, r in coloc.iterrows():
        chr_num = r["chr_num"]
        hit1 = str(r["hit1"])
        hit2 = str(r["hit2"])
        nfr_peak = str(r["NFR_peak"])
        nuc_peak = str(r["Nucleosomal_peak"])

        h1_nfr = nfr_hits.get((chr_num, hit1, nfr_peak), {})
        h1_nuc = nuc_hits.get((chr_num, hit1, nuc_peak), {})
        h2_nfr = nfr_hits.get((chr_num, hit2, nfr_peak), {})
        h2_nuc = nuc_hits.get((chr_num, hit2, nuc_peak), {})

        a1 = ann.get(hit1, {})
        a2 = ann.get(hit2, {})

        rows.append(
            {
                "Chr": r["Chr"],
                "hit1": hit1,
                "hit2": hit2,
                "NFR_peak": nfr_peak,
                "Nucleosomal_peak": nuc_peak,
                "h1_effect_allele": a1.get("effect_allele"),
                "h1_other_allele": a1.get("other_allele"),
                "h1_eaf": a1.get("eaf"),
                "h2_effect_allele": a2.get("effect_allele"),
                "h2_other_allele": a2.get("other_allele"),
                "h2_eaf": a2.get("eaf"),
                "h1_nfr_beta": h1_nfr.get("slope"),
                "h1_nfr_se": h1_nfr.get("slope_se"),
                "h1_nfr_pval": h1_nfr.get("nom_pval"),
                "h1_nuc_beta": h1_nuc.get("slope"),
                "h1_nuc_se": h1_nuc.get("slope_se"),
                "h1_nuc_pval": h1_nuc.get("nom_pval"),
                "h2_nfr_beta": h2_nfr.get("slope"),
                "h2_nfr_se": h2_nfr.get("slope_se"),
                "h2_nfr_pval": h2_nfr.get("nom_pval"),
                "h2_nuc_beta": h2_nuc.get("slope"),
                "h2_nuc_se": h2_nuc.get("slope_se"),
                "h2_nuc_pval": h2_nuc.get("nom_pval"),
                "samplesize_nfr": args.n_exp,
                "samplesize_nuc": args.n_out,
            }
        )

    out_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
