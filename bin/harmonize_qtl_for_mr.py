#!/usr/bin/env python3

import argparse
import os
import pandas as pd


def to_float(series):
    return pd.to_numeric(series, errors='coerce')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nfr', required=True)
    ap.add_argument('--nuc', required=True)
    ap.add_argument('--n-exp', required=True, type=int)
    ap.add_argument('--n-out', required=True, type=int)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    exp = pd.read_csv(args.nfr, sep='\t', dtype=str)
    out = pd.read_csv(args.nuc, sep='\t', dtype=str)

    exp = exp.rename(columns={
        'hit1': 'SNP',
        'slope': 'beta.exposure',
        'slope_se': 'se.exposure',
        'nom_pval': 'pval.exposure',
        'effect_allele': 'effect_allele.exposure',
        'other_allele': 'other_allele.exposure',
        'eaf': 'eaf.exposure',
    })

    out = out.rename(columns={
        'hit2': 'SNP',
        'slope': 'beta.outcome',
        'slope_se': 'se.outcome',
        'nom_pval': 'pval.outcome',
        'effect_allele': 'effect_allele.outcome',
        'other_allele': 'other_allele.outcome',
        'eaf': 'eaf.outcome',
    })

    keep_exp = [
        'SNP', 'beta.exposure', 'se.exposure', 'pval.exposure',
        'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure'
    ]
    keep_out = [
        'SNP', 'beta.outcome', 'se.outcome', 'pval.outcome',
        'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome'
    ]

    exp = exp[keep_exp].drop_duplicates(subset=['SNP'])
    out = out[keep_out].drop_duplicates(subset=['SNP'])

    df = exp.merge(out, on='SNP', how='inner')

    if df.empty:
        raise ValueError('No overlapping SNPs between exposure and outcome for harmonization.')

    df['beta.exposure'] = to_float(df['beta.exposure'])
    df['se.exposure'] = to_float(df['se.exposure'])
    df['pval.exposure'] = to_float(df['pval.exposure'])
    df['eaf.exposure'] = to_float(df['eaf.exposure'])

    df['beta.outcome'] = to_float(df['beta.outcome'])
    df['se.outcome'] = to_float(df['se.outcome'])
    df['pval.outcome'] = to_float(df['pval.outcome'])
    df['eaf.outcome'] = to_float(df['eaf.outcome'])

    ee = df['effect_allele.exposure'].astype(str)
    oe = df['other_allele.exposure'].astype(str)
    eo = df['effect_allele.outcome'].astype(str)
    oo = df['other_allele.outcome'].astype(str)

    same = (ee == eo) & (oe == oo)
    swap = (ee == oo) & (oe == eo)

    df_h = df[same | swap].copy()

    swap_idx = swap[same | swap]
    df_h.loc[swap_idx, 'beta.outcome'] = -1.0 * df_h.loc[swap_idx, 'beta.outcome']

    has_eaf = df_h['eaf.outcome'].notna()
    idx = swap_idx & has_eaf
    df_h.loc[idx, 'eaf.outcome'] = 1.0 - df_h.loc[idx, 'eaf.outcome']

    df_h['samplesize.exposure'] = args.n_exp
    df_h['samplesize.outcome'] = args.n_out

    out_cols = [
        'SNP',
        'beta.exposure', 'se.exposure', 'pval.exposure', 'effect_allele.exposure', 'other_allele.exposure', 'eaf.exposure', 'samplesize.exposure',
        'beta.outcome', 'se.outcome', 'pval.outcome', 'effect_allele.outcome', 'other_allele.outcome', 'eaf.outcome', 'samplesize.outcome',
    ]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df_h[out_cols].to_csv(args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()
