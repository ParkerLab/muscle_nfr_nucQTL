configfile: "caQTL_np.yaml"

celltypes=config['celltypes']
fpc_nums=config["fpc_num"]
gpc_nums=config["gpc_num"]
ids=config["IDS"]
bam_root=config["bam_root"]
qtltools_container=config["qtltools_container"]

rule all:
    input:
        expand("../../results/results_02262026/counts_file_tidy/figures/{celltype}_ca_narrowpeak_feature_summary.png", celltype=celltypes),
        expand("../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}_ca_F{fpc}_G{gpc}/permutation_chr22.txt", celltype=celltypes, fpc=fpc_nums, gpc=gpc_nums),
        expand("../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}_ca_F{fpc}_G{gpc}/merged_QTLresults_fdrcorr.csv", celltype=celltypes, fpc=fpc_nums, gpc=gpc_nums),
        expand("../../results/results_02262026/QTL_summary/{celltype}_opt_fpc.csv", celltype=celltypes),
        expand("../../results/results_02262026/QTL_opt_results/{celltype}_{region}/nominal_chr22.txt", celltype=celltypes, region="ca"),
        expand("../../results/results_02262026/QTL_opt_results/{celltype}_ca/merged_QTLresults_fdrcorr_q0.05.bed", celltype=celltypes),
        expand("../../results/results_02262026/{celltype}_ca_susie/chr{id}/chr{id}_susie_cset95_summary.tsv", celltype=celltypes, id=ids)
rule bedtosaf:
    input:
        caPeaks="../../results/results_08042024/macs2/{celltype}/atac_macs2_{celltype}_peaks.narrowPeak",
    params:
        celltype = "{celltype}",
    output:
        saf="../../results/results_02262026/caPeaks/ca_{celltype}_narrowpeak.saf",
    shell:
        """
        awk 'OFS="\t" {{print $1"."$2"."$3, $1, $2, $3, "."}}' {input} > {output.saf}
        """

rule countmatrix:
    input:
        "../../results/results_02262026/caPeaks/ca_{celltype}_narrowpeak.saf"
    params:
        celltype="{celltype}",
    output:
        "../../results/results_02262026/counts_file/{celltype}_narrowpeak_cacounts"
    conda:
        "atac"
    shell:
        """
        path={bam_root}/{params.celltype}
        featureCounts \
            -F SAF \
            -O \
            --minOverlap 1 \
            -T 16 \
            -p --donotsort \
            -a {input} \
            --countReadPairs \
            -o {output} \
            ${{path}}/*.bam
        """

rule tidy_rawcounts:
    input:
        "../../results/results_02262026/counts_file/{celltype}_narrowpeak_cacounts"
    params:
        celltype="{celltype}",
    output:
        "../../results/results_02262026/counts_file_tidy/{celltype}_ca_counts_narrowpeak_tidy.rds"
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

rule feature_summary:
    input:
        "../../results/results_02262026/counts_file_tidy/{celltype}_ca_counts_narrowpeak_tidy.rds"
    params:
        celltype="{celltype}",
    output:
        "../../results/results_02262026/counts_file_tidy/figures/{celltype}_ca_narrowpeak_feature_summary.png",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/feature_summary.py {params.celltype} ca {input} {output} 1 2 3 4 5 10 15 20
        """

rule feature_selection:
    input:
        "../../results/results_02262026/counts_file_tidy/{celltype}_ca_counts_narrowpeak_tidy.rds"
    params:
        celltype="{celltype}",
    output:
        "../../results/results_02262026/counts_file_filtered/{celltype}_ca_counts_narrowpeak_filtered.txt",
    conda:
        "atac"
    shell:
        """
        # for nfr 3 reads and 0.05
        python ../../scripts/bin/feature_selection.py {input} 3 0.05 {output}
        """

rule normalization:
    input:
        "../../results/results_02262026/counts_file_filtered/{celltype}_ca_counts_narrowpeak_filtered.txt"
    params:
        celltype="{celltype}",
        ext="narrowpeak"
    output:
        "../../results/results_02262026/normalized_counts/{celltype}_ca_narrowpeak_normalized.rds"
    conda:
        "R423"
    shell:
        """
        mkdir -p ../../results/results_02262026/QC/{params.celltype}_ca_{params.ext}/
        mkdir -p ../../results/results_02262026/bed_peaks_filtered_{params.ext}/
        Rscript ../../scripts/bin/normalization.R \
            --counts {input} \
            --counts_norm {output} \
            --outdirQC ../../results/results_02262026/QC/{params.celltype}_ca_{params.ext}/ \
            --outdirBED ../../results/results_02262026/bed_peaks_filtered_{params.ext}/ \
            --sample_info_file {config["sample_info_file"]} \
            --celltype {params.celltype} \
            --region ca
        """

rule counts_to_bed:
    input:
        "../../results/results_02262026/normalized_counts/{celltype}_ca_narrowpeak_normalized.rds"
    output:
        "../../results/results_02262026/QTL_files/counts_bed/{celltype}_ca_narrowpeak_normalized.bed"
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/counts_to_bed.py {input} {output}
        """

rule get_covariates:
    input:
        "../../results/results_02262026/normalized_counts/{celltype}_ca_narrowpeak_normalized.rds"
    params:
        outdir="../../results/results_02262026/covariates_narrowpeak/",
        fpc="{fpc}",
        gpc="{gpc}",
        vcf_eigenvec=config["vcf_eigenvec"],
        celltype="{celltype}",
        figdir="../../results/results_02262026/heatmaps_covariates_narrowpeak/",
        sample_info=lambda wildcards: config["sample_info_template"].format(celltype=wildcards.celltype),
        region='ca',
    output:
        "../../results/results_02262026/covariates_narrowpeak/qtl_cov_simple_{celltype}_ca_FPC{fpc}_GPC{gpc}.txt"
    conda:
        "atac"
    shell:
        """
        mkdir -p {params.figdir}
        python ../../scripts/bin/get_covariates.py {params.fpc} {params.gpc} {params.outdir} {params.figdir} {params.sample_info} {params.vcf_eigenvec} {input} {params.celltype} {params.region}
        """

rule preparebed:
    input:
        bed="../../results/results_02262026/QTL_files/counts_bed/{celltype}_ca_narrowpeak_normalized.bed"
    params:
        celltype="{celltype}",
        region = "ca",
        localbed="{celltype}_ca_narrowpeak_normalized.bed",
    output:
        "{celltype}_ca_narrowpeak_normalized.bed.gz"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        scp {input.bed} {params.localbed}
        bgzip {params.localbed} && tabix -p bed {params.localbed}.gz
        """

rule QTL_scan_nocov:
    input:
        vcf=ancient(config["vcf_file"]),
        bed="../../results/results_02262026/QTL_files/counts_bed/{celltype}_ca_narrowpeak_normalized.bed",
        cov="../../results/results_02262026/covariates_narrowpeak/qtl_cov_simple_{celltype}_ca_FPC{fpc}_GPC{gpc}.txt",
    params:
        window_size=config["window_size"],
        fpc="{fpc}",
        gpc="{gpc}",
        celltype="{celltype}",
        region = "ca",
        outdir="../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}_ca_F{fpc}_G{gpc}/",
        localbed="{celltype}_ca_F{fpc}_G{gpc}_narrowpeak_normalized.bed",
        localoutdir="{celltype}_ca_F{fpc}_G{gpc}",
        localcov="qtl_cov_simple_{celltype}_ca_FPC{fpc}_GPC{gpc}.txt",
    output:
        "../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}_ca_F{fpc}_G{gpc}/permutation_chr22.txt"
    envmodules:
        "singularity",
        "Bioinformatics",
        "samtools/1.13-fwwss5n",
    shell:
        """
        scp {input.bed} {params.localbed}
        bgzip -f {params.localbed}
        tabix -p bed {params.localbed}.gz
        scp {input.cov} .
        mkdir -p {params.localoutdir}
        mkdir -p {params.outdir}
        for i in $(seq 1 22)
        do
            singularity exec {qtltools_container} QTLtools cis \
                    --vcf {input.vcf} \
                    --bed {params.localbed}.gz \
                    --cov {params.localcov} \
                    --region chr${{i}} \
                    --permute 1000 \
                    --out {params.localoutdir}/permutation_chr${{i}}.txt \
                    --normal \
                    --window {params.window_size} \
                    --std-err \
                    --seed 123
        done
        mv {params.localoutdir}/permutation_chr*.txt {params.outdir}
        rmdir {params.localoutdir}
        rm -r {params.localcov}
        """

rule tidy_QTLresults:
    input:
        nocov="../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}_ca_F{fpc}_G{gpc}/permutation_chr22.txt",
    params:
        nocov_permutations="../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}_ca_F{fpc}_G{gpc}",
    output:
        nocov="../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}_ca_F{fpc}_G{gpc}/merged_QTLresults.txt",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/Tidy_QTLresults.py {params.nocov_permutations} 0.05
        """


rule calc_qvalue:
    input:
        nocov="../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}_ca_F{fpc}_G{gpc}/merged_QTLresults.txt",
    output:
        "../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}_ca_F{fpc}_G{gpc}/merged_QTLresults_fdrcorr.csv",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/qvalue.R {input.nocov}

        """

rule permutation_summary:
    input:
        nocov="../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}_ca_F100_G5/merged_QTLresults_fdrcorr.csv"
    output:
        "../../results/results_02262026/QTL_summary/{celltype}_qtl_by_fpcnum_summary.txt"
    params:
        outdir="../../results/results_02262026/QTL_summary/",
        celltype="{celltype}",
        nocov_dir="../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/",
        gpc="5",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/qtl_byfpc_summary_np.py {params.celltype} {params.nocov_dir} {params.outdir} {params.gpc}
        """

# preapre to run caQTL fine mapping
rule select_fpc:
    input:
        "../../results/results_02262026/QTL_summary/{celltype}_qtl_by_fpcnum_summary.txt"
    output:
        "../../results/results_02262026/QTL_summary/{celltype}_opt_fpc.csv"
    params:
        celltype="{celltype}",
    conda:
        "atac"
    shell:
        """
        # args: 1. summary txt 2. celltype 3. output
        python ../../scripts/bin/find_max_caqtl.py {input} {params.celltype} {output}
        """
rule nominal_pass:
    input:
        summary="../../results/results_02262026/QTL_summary/{celltype}_opt_fpc.csv",
        vcf=ancient(config["vcf_file"]),
        bed="../../results/results_02262026/QTL_files/counts_bed/{celltype}_ca_narrowpeak_normalized.bed",
    output:
        nominal=expand("../../results/results_02262026/QTL_opt_results/{{celltype}}_ca/nominal_chr{id}.txt", id=ids),
        cov="../../results/results_02262026/QTL_opt_results/{celltype}_ca/{celltype}_ca_simple_cov.txt",
        qtl="../../results/results_02262026/QTL_opt_results/{celltype}_ca/merged_QTLresults_fdrcorr.csv"
    params:
        celltype="{celltype}",
        region="ca",
        nocov_dir="../../results/results_02262026/QTLscan_nocov_permutation_narrowpeak/{celltype}",
        cov="../../results/results_02262026/covariates_narrowpeak/qtl_cov_simple_{celltype}_ca_",
        outdir="../../results/results_02262026/QTL_opt_results/{celltype}_ca/",
        local_nominal_dir="{celltype}_ca_nominal_tmp",
        window_size='500000',
        localcov="{celltype}_ca_simple_cov.txt",
        heatmap_dir="../../results/results_02262026/heatmaps_covariates_narrowpeak",
        local_bed="{celltype}_ca_narrowpeak_normalized.bed",
        vcf=config["vcf_file"]
    envmodules:
        "singularity",
        "Bioinformatics",
        "samtools/1.13-fwwss5n",
        "python/3.9.12"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {params.local_nominal_dir}
        scp {input.vcf}* .
        scp {input.bed} {params.local_bed}
        bgzip -f {params.local_bed}
        tabix -p bed {params.local_bed}.gz
        fpc=$(awk -F',' 'NR>1 && $1=="{params.celltype}" {{print $2; exit}}' {input.summary})
        echo ${{fpc}}
        scp {params.nocov_dir}_{params.region}_F${{fpc}}_G5/merged_QTLresults_fdrcorr.csv {output.qtl}
        # heatmaps
        scp {params.heatmap_dir}/FPC${{fpc}}_{params.celltype}_{params.region}.png {params.outdir}FPC${{fpc}}_{params.celltype}_{params.region}.png
        # get cov
        scp {params.cov}FPC${{fpc}}_GPC5.txt {output.cov}
        # local cov
        scp {params.cov}FPC${{fpc}}_GPC5.txt {params.localcov}

        for i in $(seq 1 22)
        do
            singularity exec {qtltools_container} QTLtools cis \
                    --vcf {params.vcf} \
                    --bed {params.local_bed}.gz \
                    --cov {params.localcov} \
                    --region chr${{i}} \
                    --nominal 1.0 \
                    --out {params.local_nominal_dir}/nominal_chr${{i}}.txt \
                    --normal \
                    --window {params.window_size} \
                    --std-err \
                    --seed 123
        done
        mv {params.local_nominal_dir}/nominal_chr*.txt {params.outdir}
        rmdir {params.local_nominal_dir}
        rm -r {params.localcov}*
        rm -r {params.local_bed}*
        """
rule correct_bed:
   input: 
       bed="../../results/results_02262026/QTL_files/counts_bed/{celltype}_ca_narrowpeak_normalized.bed",
       cov="../../results/results_02262026/QTL_opt_results/{celltype}_ca/{celltype}_ca_simple_cov.txt"
   params:
       localcov="{celltype}_ca_simple_cov.txt",
       localbed="{celltype}_ca_narrowpeak_normalized.bed",
       localout="{celltype}_ca_corrected.bed.gz"
   envmodules:
         "singularity",
         "Bioinformatics",
         "samtools/1.13-fwwss5n",
   output:
       "../../results/results_02262026/corrected_bed/{celltype}_ca_corrected.bed.gz"
   shell:
       """
       scp {input.bed} {params.localbed}
       bgzip -f {params.localbed}
       tabix -p bed {params.localbed}.gz
       scp {input.cov} {params.localcov}
       singularity exec {qtltools_container} QTLtools correct --bed {params.localbed}.gz \
                        --cov {params.localcov}\
                        --normal \
                        --out {params.localout}
        mv {params.localout} {output}
        rm -r {params.localcov}
        rm -r {params.localbed}*
       """

rule qtl_to_bed:
    input:
        qtl="../../results/results_02262026/QTL_opt_results/{celltype}_ca/merged_QTLresults_fdrcorr.csv"
    output:
        "../../results/results_02262026/QTL_opt_results/{celltype}_ca/merged_QTLresults_fdrcorr_q0.05.bed"
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/qtl_csv_to_bed.R \
            --csv {input.qtl} \
            --out {output}
        """

rule runsusie:
    input:
        nominal="../../results/results_02262026/QTL_opt_results/{celltype}_ca/nominal_chr{id}.txt",
        bed="../../results/results_02262026/corrected_bed/{celltype}_ca_corrected.bed.gz",
        qtl="../../results/results_02262026/QTL_opt_results/{celltype}_ca/merged_QTLresults_fdrcorr.csv"
    params:
        chr="chr{id}",
        min_corr=config["min_corr"],
        num_L=config["num_L"],
        outdir=config["outdir"],
        nominal_pass_dir=config["nominal_dir"],
        celltype="{celltype}",
        region="ca",
    output:
        "../../results/results_02262026/{celltype}_ca_susie/chr{id}/chr{id}_susie_cset95_summary.tsv"
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        mkdir -p {params.outdir}{params.celltype}_{params.region}_susie/{params.chr}/
        Rscript ../../scripts/bin/runsusie_qtl_ind.R \
            --chr {params.chr} \
            --bed {input.bed} \
            --qtl {input.qtl} \
            --min_corr {params.min_corr} \
            --num_L {params.num_L} \
            --outdir {params.outdir}{params.celltype}_{params.region}_susie/{params.chr}/ \
            --nominal {input.nominal} \
            --sample_info_file {config["sample_info_file"]} \
            --vcf_file {config["vcf_file"]} \
            --temp_vcf_dir {config["temp_vcf_dir"]}
        """
