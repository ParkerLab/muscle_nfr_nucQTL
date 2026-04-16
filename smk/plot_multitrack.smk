configfile: "plot_multitrack.yaml"

import csv
import os
import shlex


REPO_ROOT = os.path.abspath(os.path.join(workflow.basedir, "..", ".."))
BIN_DIR = os.path.join(REPO_ROOT, "scripts", "bin")


def abs_repo_path(path):
    path = str(path)
    return path if os.path.isabs(path) else os.path.join(REPO_ROOT, path)


SUMMARY_FILE = abs_repo_path(config["summary_file"])
COLOC_ROOT = abs_repo_path(config["paths"]["coloc_root"])
NFR_SUSIE_ROOT = abs_repo_path(config["paths"]["nfr_susie_root"])
NUC_SUSIE_ROOT = abs_repo_path(config["paths"]["nuc_susie_root"])
CA_SUSIE_ROOT = abs_repo_path(config["paths"]["ca_susie_root"])
GWAS_SUSIE_ROOT = abs_repo_path(config["paths"]["gwas_susie_root"])
OUTPUT_ROOT = abs_repo_path(config["paths"]["output_root"])

CELLTYPE = config["plot"]["celltype"]
GWAS_LABEL = config["plot"]["gwas_label"]
GWAS_TRAIT = config["plot"].get("gwas_trait", GWAS_LABEL)
FILTER_WINDOW = str(config["plot"]["filter_window"])

PPH4_NFR_MIN = float(config["filters"]["pph4_nfr_min"])
PPH4_NUC_MIN = float(config["filters"]["pph4_nuc_min"])
PPH4_CA_NP_MAX = float(config["filters"]["pph4_ca_np_max"])


def load_summary(path):
    def parse_value(value):
        if value in {"NA", "NaN", "nan", ""}:
            return None
        try:
            return float(value)
        except ValueError:
            return value

    with open(path) as handle:
        lines = [line.strip() for line in handle if line.strip()]

    headers = shlex.split(lines[0])
    rows = []
    for line in lines[1:]:
        values = shlex.split(line)
        if len(values) == len(headers) + 1:
            values = values[1:]
        if len(values) != len(headers):
            raise ValueError(f"Unexpected column count in {path}: {line}")
        row = {header: parse_value(value) for header, value in zip(headers, values)}
        rows.append(row)

    return rows


def keep_row(row):
    nfr = row.get("PP_H4_abf_nfr")
    nuc = row.get("PP_H4_abf_nuc")
    ca_np = row.get("PP_H4_abf_ca_np")

    passes_nfr = nfr is not None and nfr > PPH4_NFR_MIN
    passes_nuc = nuc is not None and nuc > PPH4_NUC_MIN
    passes_ca_np = ca_np is None or ca_np < PPH4_CA_NP_MAX
    return (passes_nfr or passes_nuc) and passes_ca_np


def selected_gwas_peaks():
    rows = load_summary(SUMMARY_FILE)
    peaks = sorted({row["gwas_peak"] for row in rows if keep_row(row) and row.get("gwas_peak")})
    if not peaks:
        raise ValueError("No gwas_peak values matched the configured filters.")
    return peaks


GWAS_PEAKS = selected_gwas_peaks()


rule all:
    input:
        expand(f"{OUTPUT_ROOT}/png/{{gwas_peak}}.png", gwas_peak=GWAS_PEAKS),
        f"{OUTPUT_ROOT}/selected_gwas_peaks.tsv",


rule write_peak_list:
    input:
        SUMMARY_FILE
    output:
        f"{OUTPUT_ROOT}/selected_gwas_peaks.tsv",
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        rows = [row for row in load_summary(input[0]) if keep_row(row)]
        dedup = {}
        for row in rows:
            dedup[row["gwas_peak"]] = {
                "gwas_peak": row["gwas_peak"],
                "PP_H4_abf_nfr": row.get("PP_H4_abf_nfr"),
                "PP_H4_abf_nuc": row.get("PP_H4_abf_nuc"),
                "PP_H4_abf_ca_np": row.get("PP_H4_abf_ca_np"),
            }
        selected = [dedup[key] for key in sorted(dedup)]
        with open(output[0], "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["gwas_peak", "PP_H4_abf_nfr", "PP_H4_abf_nuc", "PP_H4_abf_ca_np"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(selected)


rule plot_locuszoom:
    input:
        coloc_nfr=lambda wc: f"{COLOC_ROOT}/{CELLTYPE}_nfr/gwas_{wc.gwas_peak}_res_summary_{FILTER_WINDOW}win.Rda",
        coloc_nuc=lambda wc: f"{COLOC_ROOT}/{CELLTYPE}_nuc/gwas_{wc.gwas_peak}_res_summary_{FILTER_WINDOW}win.Rda",
        coloc_ca_np=lambda wc: f"{COLOC_ROOT}/{CELLTYPE}_ca_np/gwas_{wc.gwas_peak}_res_summary_{FILTER_WINDOW}win.Rda",
        selected=f"{OUTPUT_ROOT}/selected_gwas_peaks.tsv",
    output:
        f"{OUTPUT_ROOT}/png/{{gwas_peak}}.png",
    params:
        celltype=CELLTYPE,
        gwas_label=GWAS_LABEL,
        gwas_trait=GWAS_TRAIT,
        nfr_susie_root=NFR_SUSIE_ROOT,
        nuc_susie_root=NUC_SUSIE_ROOT,
        ca_susie_root=CA_SUSIE_ROOT,
        gwas_susie_root=GWAS_SUSIE_ROOT,
    conda:
        "R423"
    shell:
        """
        mkdir -p $(dirname {output:q})
        Rscript {BIN_DIR}/plot_multitrack_locuszoom.R \
          --gwas_peak {wildcards.gwas_peak:q} \
          --celltype {params.celltype:q} \
          --gwas_label {params.gwas_label:q} \
          --gwas_trait {params.gwas_trait:q} \
          --coloc_nfr {input.coloc_nfr:q} \
          --coloc_nuc {input.coloc_nuc:q} \
          --coloc_ca_np {input.coloc_ca_np:q} \
          --nfr_susie_root {params.nfr_susie_root:q} \
          --nuc_susie_root {params.nuc_susie_root:q} \
          --ca_susie_root {params.ca_susie_root:q} \
          --gwas_susie_root {params.gwas_susie_root:q} \
          --out_png {output:q}
        """
