configfile: "cit_match.yaml"

import os
import pandas as pd

bam_root = config.get("bam_root", "../../data/bams-by-cluster-sample")
out_root = "../../results/results_02052026/coverage_match"
downsample_root = "../../results/results_02052026/downsample_nfr_match"
counts_root = "../../results/results_02052026"

nfr_pattern = os.path.join(bam_root, "{celltype}", "nfr", "{sample}_nfr.bam")
nuc_pattern = os.path.join(bam_root, "{celltype}", "nuc", "{sample}_nuc.bam")

nfr_wc = glob_wildcards(nfr_pattern)
nuc_wc = glob_wildcards(nuc_pattern)

PPH4_threshold=config["PPH4_threshold"]
celltypes=config["celltypes"]
cit_outdir=config["outdir_cit"]
fdr_outdir=config["fdr_outdir"]
celltype_set = set(celltypes)
nfr_set = {(ct, s) for ct, s in zip(nfr_wc.celltype, nfr_wc.sample) if ct in celltype_set}
nuc_set = {(ct, s) for ct, s in zip(nuc_wc.celltype, nuc_wc.sample) if ct in celltype_set}
sample_pairs = sorted(nfr_set & nuc_set)
samples_by_celltype = {celltype: [] for celltype in celltypes}
for celltype, sample in sample_pairs:
    samples_by_celltype[celltype].append(sample)
active_celltypes = [ct for ct in celltypes if len(samples_by_celltype.get(ct, [])) > 0]

nfr_count_targets = [
    os.path.join(out_root, ct, f"{s}.nfr.readcount.txt")
    for ct, s in sample_pairs
]
nuc_count_targets = [
    os.path.join(out_root, ct, f"{s}.nuc.readcount.txt")
    for ct, s in sample_pairs
]
summary_table_by_celltype_targets = [
    os.path.join(out_root, ct, "coverage_by_sample.tsv")
    for ct in active_celltypes
]
summary_table = os.path.join(out_root, "coverage_by_sample.tsv")
downsample_bam_targets = [
    os.path.join(downsample_root, ct, f"{s}.bam")
    for ct, s in sample_pairs
]
downsample_bai_targets = [f"{p}.bai" for p in downsample_bam_targets]

def get_pairs(trait, PPH4_threshold):
    trait_summary_path = "../../results/results_07172024/coloc_plots_"
    trait_summary_path += trait + "/coloc_summary_" + trait + "_5000.tsv"
    trait_summary_df = pd.read_csv(trait_summary_path, sep = '\t', header=0)
    trait_summary_df = trait_summary_df[trait_summary_df['PP.H4.abf'] > PPH4_threshold]
    cit_pairs = []
    for i in range(len(trait_summary_df)):
        hit1 = str(trait_summary_df.iloc[i,]['hit1'])
        hit2 = str(trait_summary_df.iloc[i,]['hit2'])
        NFR_peak = str(trait_summary_df.iloc[i,]['NFR_peak'])
        Nuc_peak = str(trait_summary_df.iloc[i,]['Nucleosomal_peak'])
        cit_pair = trait + '-' + hit1 + '-' + hit2 + '-'
        cit_pair += NFR_peak + '-' + Nuc_peak
        cit_pairs.append(cit_pair)
    return(cit_pairs)

all_pairs = {celltype: get_pairs(celltype, PPH4_threshold) for celltype in active_celltypes}
cit_pairs = [pair for pairs in all_pairs.values() for pair in pairs]


def pair_celltype(cit_pair):
    for ct in active_celltypes:
        if cit_pair.startswith(f"{ct}-"):
            return ct
    raise ValueError(f"Could not infer celltype from cit_pair: {cit_pair}")


rule all:
    input:
        *nfr_count_targets,
        *nuc_count_targets,
        *summary_table_by_celltype_targets,
        summary_table,
        *downsample_bam_targets,
        *downsample_bai_targets,
        expand(f"{counts_root}/counts_file/{{celltype}}_nfr_counts", celltype=active_celltypes),
        expand(f"{counts_root}/counts_file_tidy/{{celltype}}_nfr_counts_tidy.rds", celltype=active_celltypes),
        expand(f"{counts_root}/counts_file_filtered/{{celltype}}_nfr_counts_filtered.txt", celltype=active_celltypes),
        expand(f"{counts_root}/normalized_counts/{{celltype}}_nfr_normalized.rds", celltype=active_celltypes),
        expand(f"{cit_outdir}cit_PPH4_threshold{PPH4_threshold}/{{cit_pair}}_cit_perm_cov.Rda", cit_pair=cit_pairs),
        expand(f"{fdr_outdir}{{celltype}}_citfdr_cov_{PPH4_threshold}_H4_summary.tsv", celltype=active_celltypes),


rule count_nfr_reads:
    input:
        bam=nfr_pattern
    output:
        readcount=os.path.join(out_root, "{celltype}", "{sample}.nfr.readcount.txt")
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        mkdir -p $(dirname {output.readcount})
        samtools view -c {input.bam} > {output.readcount}
        """


rule count_nuc_reads:
    input:
        bam=nuc_pattern
    output:
        readcount=os.path.join(out_root, "{celltype}", "{sample}.nuc.readcount.txt")
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        mkdir -p $(dirname {output.readcount})
        samtools view -c {input.bam} > {output.readcount}
        """


rule tabulate_read_counts:
    input:
        nfr=lambda wildcards: [
            os.path.join(out_root, wildcards.celltype, f"{s}.nfr.readcount.txt")
            for s in samples_by_celltype[wildcards.celltype]
        ],
        nuc=lambda wildcards: [
            os.path.join(out_root, wildcards.celltype, f"{s}.nuc.readcount.txt")
            for s in samples_by_celltype[wildcards.celltype]
        ]
    output:
        os.path.join(out_root, "{celltype}", "coverage_by_sample.tsv")
    run:
        def sample_from_count(path, suffix):
            name = os.path.basename(path)
            sample = name[:-len(suffix)]
            return sample

        nfr_counts = {}
        for path in input.nfr:
            key = sample_from_count(path, ".nfr.readcount.txt")
            with open(path, "r") as fh:
                nfr_counts[key] = int(fh.read().strip())

        nuc_counts = {}
        for path in input.nuc:
            key = sample_from_count(path, ".nuc.readcount.txt")
            with open(path, "r") as fh:
                nuc_counts[key] = int(fh.read().strip())

        rows = []
        for sample in samples_by_celltype[wildcards.celltype]:
            nfr_totalreads = nfr_counts[sample]
            nuc_totalreads = nuc_counts[sample]
            ratio = round((nuc_totalreads / nfr_totalreads), 2) if nfr_totalreads else float("nan")
            rows.append(
                {
                    "celltype": wildcards.celltype,
                    "sample": sample,
                    "nfr_totalreads": nfr_totalreads,
                    "nuc_totalreads": nuc_totalreads,
                    "nuc_totalreads/nfr_totalreads": ratio,
                }
            )

        df = pd.DataFrame(rows)
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        df.to_csv(output[0], sep="\t", index=False)


rule combine_read_count_tables:
    input:
        summary_table_by_celltype_targets
    output:
        summary_table
    run:
        dfs = []
        for path in input:
            if os.path.exists(path):
                dfs.append(pd.read_csv(path, sep="\t"))
        outdf = pd.concat(dfs, ignore_index=True) if len(dfs) > 0 else pd.DataFrame()
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        outdf.to_csv(output[0], sep="\t", index=False)


rule downsample_nfr_by_match_ratio:
    input:
        bam=nfr_pattern,
        table=os.path.join(out_root, "{celltype}", "coverage_by_sample.tsv")
    output:
        bam=os.path.join(downsample_root, "{celltype}", "{sample}.bam"),
        bai=os.path.join(downsample_root, "{celltype}", "{sample}.bam.bai")
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        mkdir -p $(dirname {output.bam})
        ratio=$(awk -F '\t' -v ct="{wildcards.celltype}" -v smp="{wildcards.sample}" '
            NR==1 {{
                for (i=1; i<=NF; i++) {{
                    if ($i=="celltype") c=i;
                    if ($i=="sample") s=i;
                    if ($i=="nuc_totalreads/nfr_totalreads") r=i;
                }}
                next
            }}
            c && s && r && $c==ct && $s==smp {{ print $r; exit }}
            (!c) && s && r && $s==smp {{ print $r; exit }}
        ' {input.table})

        if [ -z "$ratio" ]; then
            echo "No ratio found for {wildcards.celltype}/{wildcards.sample} in {input.table}" >&2
            exit 1
        fi

        seed=$(python -c 'import sys; x=round(float(sys.argv[1]), 2); x=min(max(x, 0.0), 0.99); print(f"42.{{int(round(x*100)):02d}}")' "$ratio")
        samtools view -b -s "$seed" {input.bam} > {output.bam}
        samtools index {output.bam}
        """


rule countmatrix:
    input:
        saf="../../results/results_05082024/nfrPeaks/nfr_{celltype}.saf"
    params:
        celltype="{celltype}",
        region="nfr",
        output_dir=f"{counts_root}/counts_file",
        bam_list=lambda wildcards: " ".join(
            os.path.join(downsample_root, wildcards.celltype, f"{sample}.bam")
            for sample in samples_by_celltype.get(wildcards.celltype, []))
    output:
        f"{counts_root}/counts_file/{{celltype}}_{{region}}_counts"
    conda:
        "atac"
    shell:
        """
        mkdir -p {params.output_dir}
        featureCounts \
            -F SAF \
            -O \
            --minOverlap 1 \
            -T 16 \
            -p --donotsort \
            -a {input.saf} \
            --countReadPairs \
            -o {output} \
            {params.bam_list}
        """


rule tidy_rawcounts:
    input:
        f"{counts_root}/counts_file/{{celltype}}_nfr_counts"
    params:
        celltype="{celltype}",
    output:
        f"{counts_root}/counts_file_tidy/{{celltype}}_nfr_counts_tidy.rds"
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/tidy_raw_counts_downsample.R \
            --celltype {params.celltype} \
            --region nfr \
            --input_count {input} \
            --output {output}
        """


rule feature_selection:
    input:
        f"{counts_root}/counts_file_tidy/{{celltype}}_nfr_counts_tidy.rds"
    params:
        celltype="{celltype}",
        peaks="../../results/results_05082024/counts_file_filtered/{celltype}_nfr_counts_filtered.txt",
    output:
        f"{counts_root}/counts_file_filtered/{{celltype}}_nfr_counts_filtered.txt",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/feature_selection_full.py {input} {params.peaks} {output}
        """


rule normalization:
    input:
        f"{counts_root}/counts_file_filtered/{{celltype}}_nfr_counts_filtered.txt"
    params:
        celltype="{celltype}",
        region = "nfr",
    output:
        f"{counts_root}/normalized_counts/{{celltype}}_nfr_normalized.rds"
    conda:
        "R423"
    shell:
        """
        mkdir -p {counts_root}/QC/{params.celltype}_{params.region}/
        mkdir -p {counts_root}/bed_peaks_filtered/
        Rscript ../../scripts/bin/normalize_match.R \
            --counts {input} \
            --counts_norm {output} \
            --outdirQC {counts_root}/QC/{params.celltype}_{params.region}/ \
            --outdirBED {counts_root}/bed_peaks_filtered/ \
            --sample_info_file {config["normalization_sample_info"]} \
            --celltype {params.celltype} \
            --region {params.region}
        """


rule cit:
    input:
        coloc=lambda wildcards: (
            f"../../results/results_07172024/coloc_plots_{pair_celltype(wildcards.cit_pair)}/"
            f"coloc_summary_{pair_celltype(wildcards.cit_pair)}_5000.tsv"
        ),
        counts_nfr=lambda wildcards: (
            f"{counts_root}/normalized_counts/{pair_celltype(wildcards.cit_pair)}_nfr_normalized.rds"
        )
    output:
        f"{cit_outdir}cit_PPH4_threshold{PPH4_threshold}" + "/{cit_pair}_cit_perm_cov.Rda"
    params:
        PPH4_threshold=config["PPH4_threshold"],
        outdir = cit_outdir,
        pair = "{cit_pair}",
        cov = config["cov_dir"],
        count_mat_dir_nuc=config["count_mat_dir"],
        count_mat_dir_nfr=f"{counts_root}/normalized_counts/",
        sample_info_file=config["sample_info_file"],
        vcf_file=config["vcf_file"],
        temp_vcf_dir=config["temp_vcf_dir"],
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        Rscript ../../scripts/bin/cit_downsample.R \
            --pair {params.pair} \
            --outdir {params.outdir} \
            --PPH4_threshold {params.PPH4_threshold} \
            --cov {params.cov} \
            --count_mat_dir_nfr {params.count_mat_dir_nfr} \
            --count_mat_dir_nuc {params.count_mat_dir_nuc} \
            --sample_info_file {params.sample_info_file} \
            --vcf_file {params.vcf_file} \
            --temp_vcf_dir {params.temp_vcf_dir}
        """


rule cit_fdr:
    input:
        expand(f"{cit_outdir}cit_PPH4_threshold{PPH4_threshold}/{{cit_pair}}_cit_perm_cov.Rda", cit_pair=cit_pairs)
    output:
        f"{fdr_outdir}" + "{celltype}_citfdr_cov_" + f"{PPH4_threshold}" + "_H4_summary.tsv"
    params:
        PPH4_threshold=config["PPH4_threshold"],
        cit_outdir=cit_outdir,
        celltype="{celltype}",
        fdr_outdir=fdr_outdir,
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        Rscript ../../scripts/bin/cit_fdr.R \
            --celltype {params.celltype} \
            --cit_outdir {params.cit_outdir} \
            --PPH4_threshold {params.PPH4_threshold} \
            --fdr_outdir {params.fdr_outdir}
        """
