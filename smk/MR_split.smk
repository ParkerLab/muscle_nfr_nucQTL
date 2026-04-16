configfile: "MR_split.yaml"

import os

# Required config keys:
# raw_counts_nfr, raw_counts_nuc, sample_info, keep_peaks_file, celltype, outdir_root
raw_counts_nfr = config["raw_counts_nfr"]
raw_counts_nuc = config["raw_counts_nuc"]
sample_info = config["sample_info"]
keep_peaks_file = config["keep_peaks_file"]
celltype = config["celltype"]
outdir_root = config["outdir_root"]

nfr_peak_col = config.get("nfr_peak_col", "NFR_peak")
nuc_peak_col = config.get("nuc_peak_col", "Nucleosomal_peak")
subsetA_batches = config.get("subsetA_batches", ["499-NM", "530-NM", "564-NM", "625-NM", "678-NM"])
exclude_donors = config.get("exclude_donors", ["donor1", "donor2", "donor3"])

normalization_script = config.get(
    "normalization_script",
    "bin/normalization.R",
)
split_helper = config.get(
    "split_helper",
    "bin/split_counts.py",
)
subset_norm_helper = config.get(
    "subset_norm_helper",
    "bin/subset_norm_by_peaks.R",
)
counts_to_bed_script = config.get(
    "counts_to_bed_script",
    "bin/counts_to_bed.py",
)
get_covariates_script = config.get(
    "get_covariates_script",
    "bin/get_covariates_MR.py",
)
extract_nominal_hits_script = config.get(
    "extract_nominal_hits_script",
    "bin/extract_coloc_nominal_hits.py",
)
annotate_vcf_script = config.get(
    "annotate_vcf_script",
    "bin/annotate_qtl_with_vcf.py",
)
harmonize_mr_script = config.get(
    "harmonize_mr_script",
    "bin/harmonize_qtl_for_mr.py",
)
harmonize_twosamplemr_script = config.get(
    "harmonize_twosamplemr_script",
    "bin/harmonize_twosamplemr.R",
)
steiger_script = config.get(
    "steiger_script",
    "bin/run_steiger_from_harmonized.R",
)
prepare_rowwise_steiger_script = config.get(
    "prepare_rowwise_steiger_script",
    "bin/prepare_rowwise_steiger_inputs.py",
)
run_rowwise_steiger_script = config.get(
    "run_rowwise_steiger_script",
    "bin/run_rowwise_steiger_directionality.R",
)
plot_rowwise_steiger_script = config.get(
    "plot_rowwise_steiger_script",
    "bin/plot_rowwise_steiger_directions.R",
)
vcf_file = config["vcf_file"]
vcf_eigenvec = config["vcf_eigenvec"]
qtltools_container = config.get("qtltools_container", "path/to/qtltools.sif")
qtl_window_size = int(config.get("qtl_window_size", 500000))
qtl_gpc_num = int(config.get("qtl_gpc_num", 5))
sample_size_exposure = int(config.get("sample_size_exposure", 140))
sample_size_outcome = int(config.get("sample_size_outcome", 141))
steiger_pval_threshold = float(config.get("steiger_pval_threshold", 5e-6))
rowwise_steiger_q_threshold = float(config.get("rowwise_steiger_q_threshold", config.get("rowwise_steiger_p_threshold", 0.05)))

split_dir = os.path.join(outdir_root, "splits", celltype)

counts_splitA = os.path.join(split_dir, "counts_splitA.tsv")
counts_splitB = os.path.join(split_dir, "counts_splitB.tsv")

qc_dir_A = os.path.join(outdir_root, "qc", celltype, "splitA_nfr")
qc_dir_B = os.path.join(outdir_root, "qc", celltype, "splitB_nuc")
bed_dir_A = os.path.join(outdir_root, "bed", celltype, "splitA_nfr")
bed_dir_B = os.path.join(outdir_root, "bed", celltype, "splitB_nuc")

rds_splitA = os.path.join(outdir_root, "rds", celltype, "splitA_nfr", "counts_norm.rds")
rds_splitB = os.path.join(outdir_root, "rds", celltype, "splitB_nuc", "counts_norm.rds")

rds_splitA_coloc = os.path.join(outdir_root, "rds", celltype, "splitA_nfr", "counts_norm_coloc.rds")
rds_splitB_coloc = os.path.join(outdir_root, "rds", celltype, "splitB_nuc", "counts_norm_coloc.rds")

qtl_dir = os.path.join(outdir_root, "qtlscan", celltype)
bed_dir_qtl = os.path.join(qtl_dir, "counts_bed")
cov_dir = os.path.join(qtl_dir, "covariates")
nominal_dir_A = os.path.join(qtl_dir, "QTL_opt_results", "splitA_nfr")
nominal_dir_B = os.path.join(qtl_dir, "QTL_opt_results", "splitB_nuc")

qtl_bed_A = os.path.join(bed_dir_qtl, "splitA_nfr_normalized.bed")
qtl_bed_B = os.path.join(bed_dir_qtl, "splitB_nuc_normalized.bed")
qtl_bed_gz_A = qtl_bed_A + ".gz"
qtl_bed_gz_B = qtl_bed_B + ".gz"
qtl_bed_tbi_A = qtl_bed_gz_A + ".tbi"
qtl_bed_tbi_B = qtl_bed_gz_B + ".tbi"

cov_simple_A = os.path.join(cov_dir, f"qtl_cov_simple_{celltype}_nfr_FPC40_GPC{qtl_gpc_num}.txt")
cov_simple_B = os.path.join(cov_dir, f"qtl_cov_simple_{celltype}_nuc_FPC20_GPC{qtl_gpc_num}.txt")

nominal_out_A = os.path.join(nominal_dir_A, "nominal_chr22.txt")
nominal_out_B = os.path.join(nominal_dir_B, "nominal_chr22.txt")
steiger_dir = os.path.join(qtl_dir, "steiger_inputs")
nfr_nominal_hits = os.path.join(steiger_dir, "nfrQTL_from_splitA_nominal.tsv")
nuc_nominal_hits = os.path.join(steiger_dir, "nucQTL_from_splitB_nominal.tsv")
nfr_annotated = os.path.join(steiger_dir, "nfrQTL_from_splitA_annotated.tsv")
nuc_annotated = os.path.join(steiger_dir, "nucQTL_from_splitB_annotated.tsv")
mr_harmonized = os.path.join(steiger_dir, "mr_harmonized.tsv")
mr_harmonized_twosamplemr = os.path.join(steiger_dir, "mr_harmonized_twosamplemr.tsv")
steiger_results = os.path.join(steiger_dir, "steiger_directionality.tsv")
rowwise_steiger_input = os.path.join(steiger_dir, "rowwise_steiger_input.tsv")
rowwise_steiger_results = os.path.join(steiger_dir, "steiger_directionality_rowwise.tsv")
steiger_plot_dir = os.path.join(qtl_dir, "figures")
rowwise_steiger_summary = os.path.join(steiger_plot_dir, "steiger_directionality_rowwise_counts.tsv")
rowwise_steiger_plot = os.path.join(steiger_plot_dir, "steiger_directionality_rowwise_barplot.png")
rowwise_steiger_hit1_summary = os.path.join(steiger_plot_dir, "steiger_directionality_rowwise_hit1_counts.tsv")
rowwise_steiger_hit1_plot = os.path.join(steiger_plot_dir, "steiger_directionality_rowwise_hit1_barplot.png")

rule all:
    input:
        counts_splitA,
        counts_splitB,
        rds_splitA,
        rds_splitB,
        rds_splitA_coloc,
        rds_splitB_coloc,
        nominal_out_A,
        nominal_out_B,
        nfr_nominal_hits,
        nuc_nominal_hits,
        nfr_annotated,
        nuc_annotated,
        mr_harmonized,
        mr_harmonized_twosamplemr,
        steiger_results,
        rowwise_steiger_input,
        rowwise_steiger_results,
        rowwise_steiger_summary,
        rowwise_steiger_plot,
        rowwise_steiger_hit1_summary,
        rowwise_steiger_hit1_plot,

# splitA is nfr, restricted to user-defined subsetA batches.
rule split_counts_nfr_to_splitA:
    input:
        raw_counts=raw_counts_nfr,
        sample_info=sample_info,
    output:
        counts_splitA,
        os.path.join(split_dir, "splitA_nfr.log"),
    params:
        subsetA_batches=subsetA_batches,
        exclude_donors=exclude_donors,
    shell:
        r"""
        mkdir -p {split_dir}
        python {split_helper} \
            --raw-counts {input.raw_counts} \
            --sample-info {input.sample_info} \
            --out {output[0]} \
            --subsetA-batches {params.subsetA_batches} \
            --subset-label A \
            --exclude-donors {params.exclude_donors} \
            > {output[1]} 2>&1
        """

# splitB is nuc, restricted to batches not in subsetA.
rule split_counts_nuc_to_splitB:
    input:
        raw_counts=raw_counts_nuc,
        sample_info=sample_info,
    output:
        counts_splitB,
        os.path.join(split_dir, "splitB_nuc.log"),
    params:
        subsetA_batches=subsetA_batches,
        exclude_donors=exclude_donors,
    conda:
        "atac"
    shell:
        r"""
        mkdir -p {split_dir}
        python {split_helper} \
            --raw-counts {input.raw_counts} \
            --sample-info {input.sample_info} \
            --out {output[0]} \
            --subsetA-batches {params.subsetA_batches} \
            --subset-label B \
            --exclude-donors {params.exclude_donors} \
            > {output[1]} 2>&1
        """

rule normalize_splitA:
    input:
        counts=counts_splitA,
        sample_info=sample_info,
    output:
        rds=rds_splitA,
        qc_done=os.path.join(qc_dir_A, "QC_done.txt"),
        bed_done=os.path.join(bed_dir_A, "BED_done.txt"),
    params:
        qc_dir=qc_dir_A,
        bed_dir=bed_dir_A,
        qc_dir_arg=qc_dir_A + "/",
        bed_dir_arg=bed_dir_A + "/",
        rds_dir=lambda wildcards, output: os.path.dirname(output.rds),
    conda:
        "R423"
    shell:
        r"""
        mkdir -p {params.qc_dir} {params.bed_dir} {params.rds_dir}
        Rscript {normalization_script} \
            --counts {input.counts} \
            --counts_norm {output.rds} \
            --outdirQC {params.qc_dir_arg} \
            --outdirBED {params.bed_dir_arg} \
            --sample_info_file {input.sample_info} \
            --celltype {celltype} \
            --region nfr
        echo done > {output.qc_done}
        echo done > {output.bed_done}
        """

rule normalize_splitB:
    input:
        counts=counts_splitB,
        sample_info=sample_info,
    output:
        rds=rds_splitB,
        qc_done=os.path.join(qc_dir_B, "QC_done.txt"),
        bed_done=os.path.join(bed_dir_B, "BED_done.txt"),
    params:
        qc_dir=qc_dir_B,
        bed_dir=bed_dir_B,
        qc_dir_arg=qc_dir_B + "/",
        bed_dir_arg=bed_dir_B + "/",
        rds_dir=lambda wildcards, output: os.path.dirname(output.rds),
    conda:
        "R423"
    shell:
        r"""
        mkdir -p {params.qc_dir} {params.bed_dir} {params.rds_dir}
        Rscript {normalization_script} \
            --counts {input.counts} \
            --counts_norm {output.rds} \
            --outdirQC {params.qc_dir_arg} \
            --outdirBED {params.bed_dir_arg} \
            --sample_info_file {input.sample_info} \
            --celltype {celltype} \
            --region nuc
        echo done > {output.qc_done}
        echo done > {output.bed_done}
        """

# Subset normalized matrices to coloc peak sets after normalization.
rule subset_norm_splitA_peaks:
    input:
        norm_rds=rds_splitA,
        keep_peaks_file=keep_peaks_file,
    output:
        rds_splitA_coloc,
    params:
        peak_col=nfr_peak_col,
    conda:
        "R423"
    shell:
        r"""
        Rscript {subset_norm_helper} \
            --norm_rds {input.norm_rds} \
            --keep_peaks_file {input.keep_peaks_file} \
            --keep_peaks_col {params.peak_col} \
            --out_rds {output}
        """

rule subset_norm_splitB_peaks:
    input:
        norm_rds=rds_splitB,
        keep_peaks_file=keep_peaks_file,
    output:
        rds_splitB_coloc,
    params:
        peak_col=nuc_peak_col,
    conda:
        "R423"
    shell:
        r"""
        Rscript {subset_norm_helper} \
            --norm_rds {input.norm_rds} \
            --keep_peaks_file {input.keep_peaks_file} \
            --keep_peaks_col {params.peak_col} \
            --out_rds {output}
        """

rule counts_to_bed_splitA:
    input:
        rds_splitA
    output:
        qtl_bed_A
    conda:
        "atac"
    shell:
        r"""
        mkdir -p {bed_dir_qtl}
        python {counts_to_bed_script} {input} {output}
        """

rule counts_to_bed_splitB:
    input:
        rds_splitB
    output:
        qtl_bed_B
    conda:
        "atac"
    shell:
        r"""
        mkdir -p {bed_dir_qtl}
        python {counts_to_bed_script} {input} {output}
        """

rule preparebed_splitA:
    input:
        qtl_bed_A
    output:
        gz=qtl_bed_gz_A,
        tbi=qtl_bed_tbi_A
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        r"""
        bgzip -f -c {input} > {output.gz}
        tabix -f -p bed {output.gz}
        """

rule preparebed_splitB:
    input:
        qtl_bed_B
    output:
        gz=qtl_bed_gz_B,
        tbi=qtl_bed_tbi_B
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        r"""
        bgzip -f -c {input} > {output.gz}
        tabix -f -p bed {output.gz}
        """

rule get_covariates_splitA:
    input:
        norm_rds=rds_splitA
    output:
        simple=cov_simple_A
    params:
        fpc=40,
        gpc=qtl_gpc_num,
        region="nfr"
    conda:
        "atac"
    shell:
        r"""
        mkdir -p {cov_dir}
        python {get_covariates_script} {params.fpc} {params.gpc} {output.simple} {vcf_eigenvec} {input.norm_rds}
        """

rule get_covariates_splitB:
    input:
        norm_rds=rds_splitB
    output:
        simple=cov_simple_B
    params:
        fpc=20,
        gpc=qtl_gpc_num,
        region="nuc"
    conda:
        "atac"
    shell:
        r"""
        mkdir -p {cov_dir}
        python {get_covariates_script} {params.fpc} {params.gpc} {output.simple} {vcf_eigenvec} {input.norm_rds}
        """

rule nominal_pass_splitA:
    input:
        vcf=vcf_file,
        bed=qtl_bed_gz_A,
        bed_tbi=qtl_bed_tbi_A,
        cov=cov_simple_A
    output:
        nominal_out_A
    params:
        outdir=nominal_dir_A,
        window_size=qtl_window_size,
        qtltools_container=qtltools_container
    envmodules:
        "singularity",
        "Bioinformatics",
        "samtools/1.13-fwwss5n",
        "python/3.9.12"
    shell:
        r"""
        mkdir -p {params.outdir}

        for i in $(seq 1 22)
        do
            singularity exec {params.qtltools_container} QTLtools cis \
                    --vcf {input.vcf} \
                    --bed {input.bed} \
                    --cov {input.cov} \
                    --region chr${{i}} \
                    --nominal 1.0 \
                    --out {params.outdir}/nominal_chr${{i}}.txt \
                    --normal \
                    --window {params.window_size} \
                    --std-err \
                    --seed 123
        done
        """

rule nominal_pass_splitB:
    input:
        vcf=vcf_file,
        bed=qtl_bed_gz_B,
        bed_tbi=qtl_bed_tbi_B,
        cov=cov_simple_B
    output:
        nominal_out_B
    params:
        outdir=nominal_dir_B,
        window_size=qtl_window_size,
        qtltools_container=qtltools_container
    envmodules:
        "singularity",
        "Bioinformatics",
        "samtools/1.13-fwwss5n",
        "python/3.9.12"
    shell:
        r"""
        mkdir -p {params.outdir}

        for i in $(seq 1 22)
        do
            singularity exec {params.qtltools_container} QTLtools cis \
                    --vcf {input.vcf} \
                    --bed {input.bed} \
                    --cov {input.cov} \
                    --region chr${{i}} \
                    --nominal 1.0 \
                    --out {params.outdir}/nominal_chr${{i}}.txt \
                    --normal \
                    --window {params.window_size} \
                    --std-err \
                    --seed 123
        done
        """

rule extract_coloc_nominal_hits:
    input:
        coloc=keep_peaks_file,
        nominalA=nominal_out_A,
        nominalB=nominal_out_B
    output:
        nfr_nominal_hits,
        nuc_nominal_hits
    params:
        nominal_dir_a=nominal_dir_A,
        nominal_dir_b=nominal_dir_B
    conda:
        "atac"
    shell:
        r"""
        python {extract_nominal_hits_script} \
            --coloc {input.coloc} \
            --nominal-dir-a {params.nominal_dir_a} \
            --nominal-dir-b {params.nominal_dir_b} \
            --out-nfr {output[0]} \
            --out-nuc {output[1]}
        """

rule annotate_nominal_with_vcf:
    input:
        nfr=nfr_nominal_hits,
        nuc=nuc_nominal_hits,
        vcf=vcf_file
    output:
        nfr_annotated,
        nuc_annotated
    conda:
        "atac"
    shell:
        r"""
        python {annotate_vcf_script} \
            --nfr {input.nfr} \
            --nuc {input.nuc} \
            --vcf {input.vcf} \
            --out-nfr {output[0]} \
            --out-nuc {output[1]}
        """

rule harmonize_for_mr:
    input:
        nfr=nfr_annotated,
        nuc=nuc_annotated
    output:
        mr_harmonized
    conda:
        "atac"
    params:
        n_exp=sample_size_exposure,
        n_out=sample_size_outcome
    shell:
        r"""
        python {harmonize_mr_script} \
            --nfr {input.nfr} \
            --nuc {input.nuc} \
            --n-exp {params.n_exp} \
            --n-out {params.n_out} \
            --out {output}
        """

rule harmonize_for_mr_twosamplemr:
    input:
        nfr=nfr_annotated,
        nuc=nuc_annotated
    output:
        mr_harmonized_twosamplemr
    conda:
        "R423"
    params:
        n_exp=sample_size_exposure,
        n_out=sample_size_outcome
    shell:
        r"""
        Rscript {harmonize_twosamplemr_script} \
            --nfr {input.nfr} \
            --nuc {input.nuc} \
            --n_exp {params.n_exp} \
            --n_out {params.n_out} \
            --out {output}
        """

rule steiger_directionality:
    input:
        mr_harmonized_twosamplemr
    output:
        steiger_results
    conda:
        "R423"
    params:
        pval_threshold=steiger_pval_threshold
    shell:
        r"""
        Rscript {steiger_script} \
            --harmonized {input} \
            --pval_threshold {params.pval_threshold} \
            --out {output}
        """


rule prepare_rowwise_steiger_inputs:
    input:
        coloc=keep_peaks_file,
        nominalA=nominal_out_A,
        nominalB=nominal_out_B,
        vcf=vcf_file
    output:
        rowwise_steiger_input
    params:
        nominal_dir_a=nominal_dir_A,
        nominal_dir_b=nominal_dir_B,
        n_exp=sample_size_exposure,
        n_out=sample_size_outcome
    conda:
        "atac"
    shell:
        r"""
        python {prepare_rowwise_steiger_script} \
            --coloc {input.coloc} \
            --nominal-dir-a {params.nominal_dir_a} \
            --nominal-dir-b {params.nominal_dir_b} \
            --vcf {input.vcf} \
            --n-exp {params.n_exp} \
            --n-out {params.n_out} \
            --out {output}
        """

rule steiger_directionality_rowwise:
    input:
        rowwise_steiger_input
    output:
        rowwise_steiger_results
    conda:
        "R423"
    params:
        q_threshold=rowwise_steiger_q_threshold
    shell:
        r"""
        Rscript {run_rowwise_steiger_script} \
            --input {input} \
            --q-threshold {params.q_threshold} \
            --out {output}
        """

rule plot_rowwise_steiger_directionality:
    input:
        rowwise_steiger_results
    output:
        png=rowwise_steiger_plot,
        summary=rowwise_steiger_summary
    params:
        celltype=celltype,
        outdir=steiger_plot_dir
    conda:
        "R423"
    shell:
        r"""
        mkdir -p {params.outdir}
        Rscript {plot_rowwise_steiger_script} \
            --input {input} \
            --celltype {params.celltype} \
            --hits both \
            --out {output.png} \
            --summary {output.summary}
        """

rule plot_rowwise_steiger_directionality_hit1:
    input:
        rowwise_steiger_results
    output:
        png=rowwise_steiger_hit1_plot,
        summary=rowwise_steiger_hit1_summary
    params:
        celltype=celltype,
        outdir=steiger_plot_dir
    conda:
        "R423"
    shell:
        r"""
        mkdir -p {params.outdir}
        Rscript {plot_rowwise_steiger_script} \
            --input {input} \
            --celltype {params.celltype} \
            --hits hit1 \
            --out {output.png} \
            --summary {output.summary}
        """
