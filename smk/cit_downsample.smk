configfile: "cit.yaml"

filter_window=config["filter_window"]
PPH4_threshold=config["PPH4_threshold"]
#traits=config["traits"]
celltypes=config['celltypes']

import os
import pandas as pd

bam_root = config.get("bam_root", "../../data/bams-by-cluster-sample")
downsample_root = config.get("downsample_root", "../../data/downsample_nfr")
downsample_tags = {"ds25": 0.25, "ds50": 0.50, "ds75": 0.75}

bam_pattern = os.path.join(bam_root, "{celltype}", "nfr", "{sample}.bam")
bam_samples = glob_wildcards(bam_pattern)
samples_by_celltype = {celltype: [] for celltype in celltypes}
for celltype, sample in zip(bam_samples.celltype, bam_samples.sample):
    if celltype in samples_by_celltype:
        samples_by_celltype[celltype].append(sample)

downsample_targets = []
for celltype, sample_list in samples_by_celltype.items():
    for sample in sample_list:
        for tag in downsample_tags:
            base = os.path.join(downsample_root, celltype, f"{sample}.{tag}")
            downsample_targets.extend([f"{base}.bam", f"{base}.bam.bai"])

def get_pairs(trait, PPH4_threshold):
    trait_summary_path = "../../results/results_07172024/coloc_plots_"
    trait_summary_path += trait + "/coloc_summary_" + trait + "_5000.tsv"
    trait_summary_df = pd.read_csv(trait_summary_path, sep = '\t', header=0)
    # filter by PPH4
    trait_summary_df = trait_summary_df[trait_summary_df['PP.H4.abf'] > PPH4_threshold]
    # init lists
    cit_pairs = []
    for i in range(len(trait_summary_df)):
        hit1 = str(trait_summary_df.iloc[i,]['hit1'])
        hit2 = str(trait_summary_df.iloc[i,]['hit2'])
        NFR_peak = str(trait_summary_df.iloc[i,]['NFR_peak'])
        Nuc_peak = str(trait_summary_df.iloc[i,]['Nucleosomal_peak'])
        # paste all information in
        cit_pair = trait + '-' + hit1 + '-' + hit2 + '-'
        cit_pair += NFR_peak + '-' + Nuc_peak
        # append to lsit
        cit_pairs.append(cit_pair)
    return(cit_pairs)

# Get the pairs for each celltype
all_pairs = {celltype: get_pairs(celltype, PPH4_threshold) for celltype in celltypes}

# Create a list of all pair files for the rule input
cit_pairs = [pair for pairs in all_pairs.values() for pair in pairs]

#cit_pairs = cit_pairs

rule all:
    input:
        #expand("../../results/results_05222024/cit_PPH4_threshold{PPH4_threshold}/{cit_pair}_cit_perm_cov.Rda", PPH4_threshold=PPH4_threshold, cit_pair=cit_pairs)
        #expand("../../results/results_07212024/{celltype}_citfdr_cov_{PPH4_threshold}_H4_summary.tsv", celltype=celltypes, PPH4_threshold=PPH4_threshold),
        *downsample_targets,
        expand("../../results/results_01222026/counts_file/{celltype}_nfr_{tag}_counts", celltype=celltypes, tag=downsample_tags.keys()),
        expand("../../results/results_01222026/counts_file_tidy_{tag}/{celltype}_nfr_counts_tidy.rds", celltype=celltypes, tag=downsample_tags.keys()),
        expand("../../results/results_01222026/counts_file_filtered_{tag}/{celltype}_nfr_counts_filtered.txt", celltype=celltypes, tag=downsample_tags.keys()),
        expand("../../results/results_01222026/normalized_counts_{tag}/{celltype}_nfr_normalized.rds", celltype=celltypes, tag=downsample_tags.keys()),
        expand("../../results/results_01262026/{tag}_cit_PPH4_threshold{PPH4_threshold}/{cit_pair}_cit_perm_cov.Rda", PPH4_threshold=PPH4_threshold, cit_pair=cit_pairs, tag=downsample_tags.keys()),
        expand("../../results/results_01262026/{tag}_{celltype}_citfdr_cov_{PPH4_threshold}_H4_summary.tsv", celltype=celltypes, PPH4_threshold=PPH4_threshold, tag=downsample_tags.keys()),

rule downsample:
    wildcard_constraints:
        tag="|".join(downsample_tags.keys())
    input:
        bam=os.path.join(bam_root, "{celltype}", "nfr", "{sample}.bam")
    output:
        bam=os.path.join(downsample_root, "{celltype}", "{sample}.{tag}.bam"),
        bai=os.path.join(downsample_root, "{celltype}", "{sample}.{tag}.bam.bai")
    params:
        downsample_root=downsample_root,
        sampling_seed=lambda wildcards: f"42.{int(downsample_tags[wildcards.tag] * 100):02d}"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        mkdir -p {params.downsample_root}/{wildcards.celltype}
        samtools view -b -s {params.sampling_seed} {input.bam} > {output.bam}
        samtools index {output.bam}
        """

rule countmatrix:
    wildcard_constraints:
        tag="|".join(downsample_tags.keys())
    input:
        saf="../../results/results_05082024/nfrPeaks/nfr_{celltype}.saf"
    params:
        celltype="{celltype}",
        region="nfr",
        output_dir="../../results/results_01222026/counts_file",
        bam_list=lambda wildcards: " ".join(
            os.path.join(downsample_root, wildcards.celltype, f"{sample}.{wildcards.tag}.bam")
            for sample in samples_by_celltype.get(wildcards.celltype, []))
    output:
        "../../results/results_01222026/counts_file/{celltype}_{region}_{tag}_counts"
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
        "../../results/results_01222026/counts_file/{celltype}_nfr_{tag}_counts"
    params:
        celltype="{celltype}",
        #region = "{region}"
    output:
        #"../../results/results_05082024/counts_file_tidy/{celltype}_nfr_counts_tidy.rds"
        "../../results/results_01222026/counts_file_tidy_{tag}/{celltype}_nfr_counts_tidy.rds"
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
        "../../results/results_01222026/counts_file_tidy_{tag}/{celltype}_nfr_counts_tidy.rds"
    params:
        celltype="{celltype}",
        #region = "{region}"
        peaks="../../results/results_05082024/counts_file_filtered/{celltype}_nfr_counts_filtered.txt",
    output:
        "../../results/results_01222026/counts_file_filtered_{tag}/{celltype}_nfr_counts_filtered.txt",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/feature_selection_full.py {input} {params.peaks} {output}
        """

rule normalization:
    input:
        "../../results/results_01222026/counts_file_filtered_{tag}/{celltype}_nfr_counts_filtered.txt"
    params:
        celltype="{celltype}",
        region = "nfr",
        tag="{tag}"

    output:
        "../../results/results_01222026/normalized_counts_{tag}/{celltype}_nfr_normalized.rds"
    conda:
        "R423"
    shell:
        """
        mkdir -p ../../results/results_01222026/QC_{tag}/{params.celltype}_{params.region}/
        mkdir -p ../../results/results_01222026/bed_peaks_filtered_{tag}/
        Rscript ../../scripts/bin/normalization.R \
            --counts {input} \
            --counts_norm {output} \
            --outdirQC ../../results/results_01222026/QC_{tag}/{params.celltype}_{params.region}/ \
            --outdirBED ../../results/results_01222026/bed_peaks_filtered_{tag}/ \
            --sample_info_file {config["normalization_sample_info"]} \
            --celltype {params.celltype} \
            --region {params.region}
        """

rule cit:
    input:
        expand("../../results/results_07172024/coloc_plots_{celltype}/coloc_summary_{celltype}_5000.tsv", celltype=celltypes)
    output:
        "../../results/results_01262026/{tag}_cit_PPH4_threshold{PPH4_threshold}/{cit_pair}_cit_perm_cov.Rda"
    params:
        PPH4_threshold=config["PPH4_threshold"],
        outdir = config["outdir_cit"],
        #trait = lambda wildcards: '{trait}',
        #celltype = "{celltype}",
        pair = "{cit_pair}",
        cov = config["cov_dir"],
        count_mat_dir_nuc=config["count_mat_dir"],
        count_mat_dir_nfr="../../results/results_01222026/normalized_counts_{tag}/",
        tag="{tag}",
        sample_info_file=config["sample_info_file"],
        vcf_file=config["vcf_file"],
        temp_vcf_dir=config["temp_vcf_dir"]
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
            --outdir {params.outdir}{params.tag}_ \
            --PPH4_threshold {params.PPH4_threshold} \
            --cov {params.cov} \
            --count_mat_dir_nfr {params.count_mat_dir_nfr} \
            --count_mat_dir_nuc {params.count_mat_dir_nuc} \
            --sample_info_file {params.sample_info_file} \
            --vcf_file {params.vcf_file} \
            --temp_vcf_dir {params.temp_vcf_dir}
        """

rule cit_fdr:
    wildcard_constraints:
        tag="|".join(downsample_tags.keys())
    input:
        expand("../../results/results_01262026/{tag}_cit_PPH4_threshold{PPH4_threshold}/{cit_pair}_cit_perm_cov.Rda", cit_pair=cit_pairs, PPH4_threshold=PPH4_threshold, tag=downsample_tags.keys())
    output:
        "../../results/results_01262026/{tag}_{celltype}_citfdr_cov_{PPH4_threshold}_H4_summary.tsv"
    params:
        PPH4_threshold=config["PPH4_threshold"],
        cit_outdir=config["outdir_cit"],
        celltype="{celltype}",
        fdr_outdir=config["fdr_outdir"],
        tag="{tag}"
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
            --cit_outdir {params.cit_outdir}{params.tag}_ \
            --PPH4_threshold {params.PPH4_threshold} \
            --fdr_outdir {params.fdr_outdir}{params.tag}_
        """
