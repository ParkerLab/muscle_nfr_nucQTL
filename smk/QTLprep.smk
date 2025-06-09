configfile: "QTLprep.yaml"

celltypes=config['celltypes']
regions=config['regions']

# get the samples
import pandas as pd

def get_samples(sample_file):
    sample_df = pd.read_csv(sample_file, sep = '\t')
    # get ids
    samples = sample_df['SNG.1ST']+'--'+sample_df['batch']
    # get unique ids
    samples = list(samples.drop_duplicates())
    # remove '-NM'
    return([s.strip('-NM') for s in samples])

samples = get_samples("../../data/sample_info/sample_level_covariates_atac.tsv")

rule all:
    input:
        #expand("../../results/results_05082024/{region}Peaks/{region}_{celltype}.saf", celltype=celltypes, region=regions),
        #expand("../../results/results_05082024/{region}Peaks/{region}_{celltype}_wname.bed", celltype=celltypes, region=regions)
        #expand("../../results/results_05282024/counts_file_byfullbam/{celltype}_nfr_counts", celltype=celltypes),
        #expand("../../results/results_05292024/counts_file_tidy_byfullbam/{celltype}_nfr_counts_tidy.rds", celltype=celltypes),
        #expand("../../results/results_05082024/counts_file_tidy/figures/nfr/{celltype}_nfr_feature_summary_byreads.png", celltype=celltypes),
        #expand("../../results/results_07122024/counts_file_tidy/figures/nuc/{celltype}_nuc_feature_summary_bymean.png", celltype=celltypes),
        #expand("../../results/results_05292024/counts_file_tidy/figures/{celltype}_nfr_feature_summary_byfullbam.png", celltype=celltypes)
        #expand("../../results/results_06042024/counts_file_filtered_byfullbam/{celltype}_nfr_counts_filtered.txt", celltype=celltypes)
        #expand("../../results/results_06042024/normalized_counts_byfullbam/{celltype}_nfr_normalized.rds", celltype=celltypes),
        #expand("../../results/results_06042024/QTL_files/counts_bed_byfullbam/{celltype}_nfr_normalized.bed", celltype=celltypes)
        #expand("../../results/results_07122024/counts_file_filtered/{celltype}_nuc_counts_filtered.txt", celltype=celltypes)
        #expand("../../results/results_07152024/normalized_counts/{celltype}_{region}_normalized.rds", celltype=celltypes, region=regions),
        expand("../../results/results_07152024/QTL_files/counts_bed/{celltype}_{region}_normalized.bed", celltype=celltypes, region=regions)
        #expand("../../results/results_04302024/counts_file_sample/{celltype}_nfr/{sample}_summary.tab", celltype=celltypes, sample=samples),
        #expand("../../results/results_07122024/counts_file_sample/{celltype}_nuc/{sample}_summary.tab", celltype=celltypes, sample=samples),
        #expand("../../results/results_07122024/counts_file_tidy/{celltype}_nuc_mean_tidy.rds", celltype=celltypes),
        #expand("../../results/results_05082024/counts_file_tidy/{celltype}_nuc_max_tidy.rds", celltype=celltypes)

rule bedtosaf:
    input:
        #"../../results/results_12192023/nucleoatac/{celltype}/{celltype}.{region}.bed"
        "../../results/results_07112024/{region}Peaks/{region}_{celltype}.bed",
    params:
        celltype = "{celltype}",
        region = "{region}"
    output:
        saf="../../results/results_07122024/{region}Peaks/{region}_{celltype}.saf",
        bed4="../../results/results_07122024/{region}Peaks/{region}_{celltype}_wname.bed"
    shell:
        """
        #awk '{{print $1"\t"$2"\t"$3}}'
        awk 'OFS="\t" {{print $1"."$2"."$3, $1, $2, $3, "."}}' {input} > {output.saf}
        awk 'OFS="\t" {{print $1, $2, $3, $1"."$2"."$3}}' {input} > {output.bed4}
        """

rule countmatrix:
    input:
        "../../results/results_07122024/{region}Peaks/{region}_{celltype}.saf"
    params:
        celltype="{celltype}",
        region = "{region}"
    output:
        "../../results/results_07122024/counts_file/{celltype}_{region}counts"
        #"../../results/results_05282024/counts_file_byfullbam/{celltype}_nfr_counts"
    conda:
        "atac"
    shell:
        """
        # path=../../data/bams-by-cluster-sample/{params.celltype}/nfr
        path=../../data/bams-by-cluster-sample/{params.celltype}
        featureCounts \
            -F SAF \
            -O \
            --minOverlap 1 \
            -T 16 \
            -p --donotsort \
            -a {input} \
            --countReadPairs \
            -o ../../results/results_07122024/counts_file/{params.celltype}_{params.region}_counts \
            ${{path}}/atac-{params.celltype}@fusion@*.bam
        """

rule countmatrix_bw:
    input:
        #nfr="../../results/results_04102024/nfrPeaks/nfr_{celltype}_wname.bed",
        nuc="../../results/results_07122024/nucPeaks/nuc_{celltype}_wname.bed",
    params:
        celltype="{celltype}",
        sample= "{sample}",
        #nfr_bw="../../data/bams-by-cluster-sample/{celltype}/nfr_bw/atac-{celltype}@fusion@{sample}_nfr.bw",
        diff_bw="../../data/bams-by-cluster-sample/{celltype}/diff_bw_mod/atac-{celltype}@fusion@{sample}_diff.bw",
    output:
        #nfr="../../results/results_04302024/counts_file_sample/{celltype}_nfr/{sample}_summary.tab",
        nuc="../../results/results_07122024/counts_file_sample/{celltype}_nuc/{sample}_summary.tab"
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n",
        "UCSC-utilities/4.0.0"
    shell:
        """
        # nucPeaks by diff.bw
        bigWigAverageOverBed {params.diff_bw} {input.nuc} {output.nuc} -minMax
        """

rule tidy_rawcounts:
    input:
        #"../../results/results_05082024/counts_file/{celltype}_nfr_counts"
        "../../results/results_05282024/counts_file_byfullbam/{celltype}_nfr_counts"
    params:
        celltype="{celltype}",
        #region = "{region}"
    output:
        #"../../results/results_05082024/counts_file_tidy/{celltype}_nfr_counts_tidy.rds"
        "../../results/results_05292024/counts_file_tidy_byfullbam/{celltype}_nfr_counts_tidy.rds"
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/tidy_raw_counts.R \
            --celltype {params.celltype} \
            --region nfr \
            --input_count {input} \
            --output {output}
        """

rule tidy_rawscore:
    input:
        "../../results/results_07122024/counts_file_sample/{celltype}_nuc/P14762--564_summary.tab"
    params:
        celltype="{celltype}",
        region = "nuc",
        indir = "../../results/results_07122024/counts_file_sample/{celltype}_nuc/",
        outdir = "../../results/results_07122024/counts_file_tidy/{celltype}_nuc_",
    output:
        rdsmax="../../results/results_07122024/counts_file_tidy/{celltype}_nuc_max_tidy.rds",
        rdsmean="../../results/results_07122024/counts_file_tidy/{celltype}_nuc_mean_tidy.rds"
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/tidy_bw_score.R \
            --celltype {params.celltype} \
            --region {params.region} \
            --input_dir {params.indir} \
            --output {params.outdir}
        """

rule feature_summary:
    input:
        #nfr="../../results/results_05082024/counts_file_tidy/{celltype}_nfr_counts_tidy.rds",
        nuc="../../results/results_07122024/counts_file_tidy/{celltype}_nuc_mean_tidy.rds",
        #full="../../results/results_05292024/counts_file_tidy_byfullbam/{celltype}_nfr_counts_tidy.rds"
    params:
        celltype="{celltype}",
        #region = "{region}"
    output:
        #nfr="../../results/results_05082024/counts_file_tidy/figures/nfr/{celltype}_nfr_feature_summary_byreads.png",
        nuc="../../results/results_07122024/counts_file_tidy/figures/nuc/{celltype}_nuc_feature_summary_bymean.png",
        #full="../../results/results_05292024/counts_file_tidy/figures/{celltype}_nfr_feature_summary_byfullbam.png"
    conda:
        "atac"
    shell:
        """
        # mkdir -p ../../results/results_04302024/counts_file_tidy/figures/
        python ../../scripts/bin/feature_summary.py {params.celltype} nuc {input.nuc} {output.nuc} 0.1 0.2 0.3 0.5 0.8 1 2 5
        """

rule feature_selection:
    input:
        #nfr="../../results/results_05082024/counts_file_tidy/{celltype}_nfr_counts_tidy.rds",
        nuc="../../results/results_07122024/counts_file_tidy/{celltype}_nuc_mean_tidy.rds"
    params:
        celltype="{celltype}",
        #region = "{region}"
    output:
        #nfr="../../results/results_05082024/counts_file_filtered/{celltype}_nfr_counts_filtered.txt",
        nuc="../../results/results_07122024/counts_file_filtered/{celltype}_nuc_counts_filtered.txt",
    conda:
        "atac"
    shell:
        """
        # for nfr 3 reads and 0.05
        python ../../scripts/bin/feature_selection.py {input.nuc} 0.5 0.05 {output.nuc}
        """

rule feature_selection_full:
    input:
        "../../results/results_05292024/counts_file_tidy_byfullbam/{celltype}_nfr_counts_tidy.rds"
    params:
        celltype="{celltype}",
        #region = "{region}"
        peaks="../../results/results_05082024/counts_file_filtered/{celltype}_nfr_counts_filtered.txt",
    output:
        "../../results/results_06042024/counts_file_filtered_byfullbam/{celltype}_nfr_counts_filtered.txt",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/feature_selection_full.py {input} {params.peaks} {output}
        """

rule normalization:
    input:
        #"../../results/results_05082024/counts_file_filtered/{celltype}_{region}_counts_filtered.txt"
        #"../../results/results_06042024/counts_file_filtered_byfullbam/{celltype}_nfr_counts_filtered.txt"
        "../../results/results_07122024/counts_file_filtered/{celltype}_nuc_counts_filtered.txt"
    params:
        celltype="{celltype}",
        region = "nuc",

    output:
        "../../results/results_07152024/normalized_counts/{celltype}_nuc_normalized.rds"
    conda:
        "R423"
    shell:
        """
        mkdir -p ../../results/results_07152024/QC/{params.celltype}_{params.region}/
        mkdir -p ../../results/results_07152024/bed_peaks_filtered/
        Rscript ../../scripts/bin/normalization.R \
            --counts {input} \
            --counts_norm {output} \
            --outdirQC ../../results/results_07152024/QC/{params.celltype}_{params.region}/ \
            --outdirBED ../../results/results_07152024/bed_peaks_filtered/ \
            --sample_info_file ../../data/sample_info/sample_level_covariates_atac.tsv \
            --celltype {params.celltype} \
            --region {params.region}
        """

rule counts_to_bed:
    input:
        #"../../results/results_05082024/normalized_counts/{celltype}_{region}_normalized.rds"
        #"../../results/results_06042024/normalized_counts_byfullbam/{celltype}_nfr_normalized.rds"
        "../../results/results_07152024/normalized_counts/{celltype}_nuc_normalized.rds"
    output:
        "../../results/results_07152024/QTL_files/counts_bed/{celltype}_nuc_normalized.bed"
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/counts_to_bed.py {input} {output}
        """