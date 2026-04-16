#!/usr/bin/env python3

import argparse
import gzip
import os
import pandas as pd


def open_maybe_gzip(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_info_field(info, key):
    prefix = key + '='
    for item in info.split(';'):
        if item.startswith(prefix):
            return item[len(prefix):]
    return None


def first_float(x):
    if x is None:
        return None
    if ',' in x:
        x = x.split(',')[0]
    try:
        return float(x)
    except Exception:
        return None


def build_targets(df_nfr, df_nuc):
    targets = set()
    for v in df_nfr['hit1'].dropna().astype(str):
        targets.add(v)
    for v in df_nuc['hit2'].dropna().astype(str):
        targets.add(v)
    return targets


def annotate_from_vcf(vcf_path, targets):
    ann = {}

    # Direct annotations for canonical CHR:POS:REF:ALT variant IDs if present.
    for v in targets:
        parts = v.split(':')
        if len(parts) == 4:
            ann[v] = {
                'variant_id': v,
                'var_chr': parts[0],
                'var_pos': parts[1],
                'other_allele': parts[2],
                'effect_allele': parts[3],
                'eaf': None,
            }

    need = {v for v in targets if v not in ann}

    with open_maybe_gzip(vcf_path) as fh:
        for line in fh:
            if not line or line[0] == '#':
                continue
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, _qual, _filt, info = cols[:8]
            alts = alt.split(',')

            af_val = parse_info_field(info, 'AF')
            af_first = first_float(af_val)

            # Match rs IDs (or any IDs in ID field)
            if vid and vid != '.':
                ids = vid.split(';')
                for one_id in ids:
                    if one_id in need:
                        ann[one_id] = {
                            'variant_id': one_id,
                            'var_chr': chrom,
                            'var_pos': pos,
                            'other_allele': ref,
                            'effect_allele': alts[0],
                            'eaf': af_first,
                        }

            # Match canonical CHR:POS:REF:ALT keys
            for i, a in enumerate(alts):
                key = f"{chrom}:{pos}:{ref}:{a}"
                if key in targets:
                    af_i = None
                    if af_val is not None:
                        parts = af_val.split(',')
                        if i < len(parts):
                            af_i = first_float(parts[i])
                    ann[key] = {
                        'variant_id': key,
                        'var_chr': chrom,
                        'var_pos': pos,
                        'other_allele': ref,
                        'effect_allele': a,
                        'eaf': af_i,
                    }

    return ann


def merge_anno(df, hit_col, ann):
    out = df.copy()
    out['variant_id'] = out[hit_col].astype(str)

    meta = pd.DataFrame.from_dict(ann, orient='index')
    if not meta.empty:
        meta = meta[['variant_id', 'var_chr', 'var_pos', 'other_allele', 'effect_allele', 'eaf']]
        out = out.merge(meta, on='variant_id', how='left')
    else:
        out['var_chr'] = None
        out['var_pos'] = None
        out['other_allele'] = None
        out['effect_allele'] = None
        out['eaf'] = None

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--nfr', required=True)
    ap.add_argument('--nuc', required=True)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--out-nfr', required=True)
    ap.add_argument('--out-nuc', required=True)
    args = ap.parse_args()

    df_nfr = pd.read_csv(args.nfr, sep='\t', dtype=str)
    df_nuc = pd.read_csv(args.nuc, sep='\t', dtype=str)

    targets = build_targets(df_nfr, df_nuc)
    ann = annotate_from_vcf(args.vcf, targets)

    nfr_annot = merge_anno(df_nfr, 'hit1', ann)
    nuc_annot = merge_anno(df_nuc, 'hit2', ann)

    os.makedirs(os.path.dirname(args.out_nfr), exist_ok=True)
    os.makedirs(os.path.dirname(args.out_nuc), exist_ok=True)

    nfr_annot.to_csv(args.out_nfr, sep='\t', index=False)
    nuc_annot.to_csv(args.out_nuc, sep='\t', index=False)


if __name__ == '__main__':
    main()
