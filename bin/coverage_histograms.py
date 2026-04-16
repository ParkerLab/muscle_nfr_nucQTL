import argparse
import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate histogram PNGs for coverage table columns."
    )
    parser.add_argument("--table", required=True, help="Path to coverage_by_sample.tsv")
    parser.add_argument("--outdir", required=True, help="Output directory for PNGs")
    parser.add_argument("--cols", required=True, help="Comma-separated list of columns")
    parser.add_argument("--bins", type=int, default=30, help="Histogram bin count")
    return parser.parse_args()


def main():
    try:
        table = snakemake.input[0]
        outdir = snakemake.params["outdir"]
        cols = snakemake.params["cols"]
        bins = 30
    except NameError:
        args = parse_args()
        table = args.table
        outdir = args.outdir
        cols = [c for c in args.cols.split(",") if c]
        bins = args.bins

    df = pd.read_csv(table, sep="\t")
    os.makedirs(outdir, exist_ok=True)

    for col in cols:
        if col not in df.columns:
            continue
        series = pd.to_numeric(df[col], errors="coerce").dropna()
        if series.empty:
            continue
        plt.figure(figsize=(6, 4))
        plt.hist(series, bins=bins, color="#3b7ddd", edgecolor="black", alpha=0.85)
        plt.title(col)
        plt.xlabel(col)
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{col}.png"), dpi=150)
        plt.close()


if __name__ == "__main__":
    main()
