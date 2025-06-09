configfile: "caQTL.yaml"

celltypes=config['celltypes']
extensions=config["extensions"]
fpc_nums=config["fpc_num"]
gpc_nums=config["gpc_num"]
ids=config["IDS"]

rule all:
    input:
        # expand("../../results/results_08292024/caPeaks/ca_{celltype}_75.bed", celltype=celltypes),
        # expand("../../results/results_08292024/caPeaks/ca_{celltype}_150.bed", celltype=celltypes),
        # expand("../../results/results_08292024/counts_file/{celltype}_{ext}_cacounts", celltype=celltypes, ext=extensions),
        # expand("../../results/results_08292024/counts_file_tidy/{celltype}_ca_counts_{ext}_tidy.rds", celltype=celltypes, ext=extensions),
        # expand("../../results/results_08302024/counts_file_tidy/figures/{celltype}_ca_{ext}_feature_summary.png", 
        #         celltype=celltypes, ext=extensions),
        # expand("../../results/results_08302024/counts_file_filtered/{celltype}_ca_counts_{ext}_filtered.txt", 
        #         celltype=celltypes, ext=extensions),
        # expand("../../results/results_08302024/QTL_files/counts_bed/{celltype}_ca_{ext}_normalized.bed",
        #         celltype=celltypes, ext=extensions),
        # expand("../../results/results_08312024/covariates_{ext}ext/qtl_cov_simple_{celltype}_ca_FPC{fpc}_GPC{gpc}.txt",
        #         celltype=celltypes, ext=extensions, fpc=fpc_nums, gpc=gpc_nums),
        # expand("../../results/results_08312024/QTLscan_nocov_permutation_{ext}ext/{celltype}_ca_F{fpc}_G{gpc}/permutation_chr22.txt",
        #         celltype=celltypes, ext=extensions, fpc=fpc_nums, gpc=gpc_nums),
        # expand("../../results/results_08312024/QTLscan_nocov_permutation_{ext}ext/{celltype}_ca_F{fpc}_G{gpc}/merged_QTLresults_fdrcorr.csv",
        #         celltype=celltypes, ext=extensions, fpc=fpc_nums, gpc=gpc_nums),
        # expand("../../results/results_08312024/QTL_summary/{celltype}_opt_fpc.csv",
        #         celltype=celltypes),
        expand("../../results/results_06032025/QTL_opt_results/{celltype}_ca/nominal_chr1.txt",
                celltype=celltypes),
        #expand("/scratch/scjp_root/scjp1/xiaoouw/Muscle_snATAC_nucQTL/results/results_06032025/corrected_bed/{celltype}_ca_corrected.bed.gz",
        #        celltype=celltypes),
        expand("../../results/results_06032025/{celltype}_ca_susie/chr{id}/chr{id}_susie_cset95_summary.tsv",
                celltype=celltypes, id=ids)

rule buildpeaks:
    input:
        ca_summits="../../results/results_08042024/macs2/{celltype}/atac_macs2_{celltype}_summits.bed",
    output:
        caPeaks_150="../../results/results_08292024/caPeaks/ca_{celltype}_150.bed",
        caPeaks_75="../../results/results_08292024/caPeaks/ca_{celltype}_75.bed",
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        # 150 ext
        awk '{{print $1"\t"$2-150"\t"$3+150}}' \
            {input.ca_summits} > {output.caPeaks_150}

        # 75 ext
        awk '{{print $1"\t"$2-75"\t"$3+75}}' \
            {input.ca_summits} > {output.caPeaks_75}
        """

rule bedtosaf:
    input:
        caPeaks="../../results/results_08292024/caPeaks/ca_{celltype}_{ext}.bed",
    params:
        celltype = "{celltype}",
    output:
        saf="../../results/results_08292024/caPeaks/ca_{celltype}_{ext}.saf",
    shell:
        """
        awk 'OFS="\t" {{print $1"."$2"."$3, $1, $2, $3, "."}}' {input} > {output.saf}
        """

rule countmatrix:
    input:
        "../../results/results_08292024/caPeaks/ca_{celltype}_{ext}.saf"
    params:
        celltype="{celltype}",
    output:
        "../../results/results_08292024/counts_file/{celltype}_{ext}_cacounts"
    conda:
        "atac"
    shell:
        """
        path=/home/xiaoouw/Muscle_snATAC_nucQTL/data/bams-by-cluster-sample/{params.celltype}
        featureCounts \
            -F SAF \
            -O \
            --minOverlap 1 \
            -T 16 \
            -p --donotsort \
            -a {input} \
            --countReadPairs \
            -o {output} \
            ${{path}}/atac-{params.celltype}@fusion@*.bam
        """

rule tidy_rawcounts:
    input:
        "../../results/results_08292024/counts_file/{celltype}_{ext}_cacounts"
    params:
        celltype="{celltype}",
    output:
        "../../results/results_08292024/counts_file_tidy/{celltype}_ca_counts_{ext}_tidy.rds"
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
        "../../results/results_08292024/counts_file_tidy/{celltype}_ca_counts_{ext}_tidy.rds"
    params:
        celltype="{celltype}",
    output:
        "../../results/results_08302024/counts_file_tidy/figures/{celltype}_ca_{ext}_feature_summary.png",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/feature_summary.py {params.celltype} ca {input} {output} 1 2 3 4 5 10 15 20
        """

rule feature_selection:
    input:
        "../../results/results_08292024/counts_file_tidy/{celltype}_ca_counts_{ext}_tidy.rds"
    params:
        celltype="{celltype}",
    output:
        "../../results/results_08302024/counts_file_filtered/{celltype}_ca_counts_{ext}_filtered.txt",
    conda:
        "atac"
    shell:
        """
        # for nfr 3 reads and 0.05
        python ../../scripts/bin/feature_selection.py {input} 3 0.05 {output}
        """

rule normalization:
    input:
        "../../results/results_08302024/counts_file_filtered/{celltype}_ca_counts_{ext}_filtered.txt"
    params:
        celltype="{celltype}",
        ext="{ext}"
    output:
        "../../results/results_08302024/normalized_counts/{celltype}_ca_{ext}_normalized.rds"
    conda:
        "R423"
    shell:
        """
        mkdir -p ../../results/results_08302024/QC/{params.celltype}_ca_{params.ext}/
        mkdir -p ../../results/results_08302024/bed_peaks_filtered_{params.ext}/
        Rscript ../../scripts/bin/normalization.R \
            --counts {input} \
            --counts_norm {output} \
            --outdirQC ../../results/results_08302024/QC/{params.celltype}_ca_{params.ext}/ \
            --outdirBED ../../results/results_08302024/bed_peaks_filtered_{params.ext}/ \
            --sample_info_file ../../data/sample_info/sample_level_covariates_atac.tsv \
            --celltype {params.celltype} \
            --region ca
        """

rule counts_to_bed:
    input:
        "../../results/results_08302024/normalized_counts/{celltype}_ca_{ext}_normalized.rds"
    output:
        "../../results/results_08302024/QTL_files/counts_bed/{celltype}_ca_{ext}_normalized.bed"
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/counts_to_bed.py {input} {output}
        """

rule get_covariates:
    input:
        "../../results/results_08302024/normalized_counts/{celltype}_ca_{ext}_normalized.rds"
    params:
        outdir="../../results/results_08312024/covariates_{ext}ext/",
        fpc="{fpc}",
        gpc="{gpc}",
        vcf_eigenvec=config["vcf_eigenvec"],
        celltype="{celltype}",
        figdir="../../results/results_08312024/heatmaps_covariates_{ext}ext/",
        sample_info="../../data/sample_info/{celltype}_sample_level_covariates_atac.tsv",
        region='ca',
    output:
        "../../results/results_08312024/covariates_{ext}ext/qtl_cov_simple_{celltype}_ca_FPC{fpc}_GPC{gpc}.txt"
    conda:
        "atac"
    shell:
        """
        mkdir -p {params.figdir}
        python ../../scripts/bin/get_covariates.py {params.fpc} {params.gpc} {params.outdir} {params.figdir} {params.sample_info} {params.vcf_eigenvec} {input} {params.celltype} {params.region}
        """

rule preparebed:
    input:
        bed="../../results/results_08302024/QTL_files/counts_bed/{celltype}_ca_{ext}_normalized.bed"
    params:
        celltype="{celltype}",
        region = "ca",
        localbed="{celltype}_ca_{ext}_normalized.bed",
    output:
        "{celltype}_ca_{ext}_normalized.bed.gz"
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
        vcf="fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz",
        bed="{celltype}_ca_{ext}_normalized.bed.gz",
        cov="../../results/results_08312024/covariates_{ext}ext/qtl_cov_simple_{celltype}_ca_FPC{fpc}_GPC{gpc}.txt",
    params:
        window_size=config["window_size"],
        fpc="{fpc}",
        gpc="{gpc}",
        celltype="{celltype}",
        region = "ca",
        outdir="../../results/results_08312024/QTLscan_nocov_permutation_{ext}ext/{celltype}_ca_F{fpc}_G{gpc}/",
        localcov="qtl_cov_simple_{celltype}_ca_FPC{fpc}_GPC{gpc}.txt",
    output:
        "../../results/results_08312024/QTLscan_nocov_permutation_{ext}ext/{celltype}_ca_F{fpc}_G{gpc}/permutation_chr22.txt"
    shell:
        """
        scp {input.cov} .
        for i in $(seq 1 22)
        do
            QTLtools cis \
                    --vcf {input.vcf} \
                    --bed {input.bed} \
                    --cov {params.localcov} \
                    --region chr${{i}} \
                    --permute 1000 \
                    --out {params.outdir}permutation_chr${{i}}.txt \
                    --normal \
                    --window {params.window_size} \
                    --std-err \
                    --seed 123
        done
        rm -r {params.localcov}
        """

rule tidy_QTLresults:
    input:
        nocov="../../results/results_08312024/QTLscan_nocov_permutation_{ext}ext/{celltype}_ca_F{fpc}_G{gpc}/permutation_chr22.txt",
    params:
        nocov_permutations="../../results/results_08312024/QTLscan_nocov_permutation_{ext}ext/{celltype}_ca_F{fpc}_G{gpc}",
    output:
        nocov="../../results/results_08312024/QTLscan_nocov_permutation_{ext}ext/{celltype}_ca_F{fpc}_G{gpc}/merged_QTLresults.txt",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/Tidy_QTLresults.py {params.nocov_permutations} 0.05
        """


rule calc_qvalue:
    input:
        nocov="../../results/results_08312024/QTLscan_nocov_permutation_{ext}ext/{celltype}_ca_F{fpc}_G{gpc}/merged_QTLresults.txt",
    output:
        "../../results/results_08312024/QTLscan_nocov_permutation_{ext}ext/{celltype}_ca_F{fpc}_G{gpc}/merged_QTLresults_fdrcorr.csv",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/qvalue.R {input.nocov}

        """
# preapre to run caQTL fine mapping
rule select_fpc:
    input:
        "../../results/results_08312024/QTL_summary/{celltype}_qtl_by_fpcnum_summary.txt"
    output:
        "../../results/results_08312024/QTL_summary/{celltype}_opt_fpc.csv"
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
        summary="../../results/results_08312024/QTL_summary/{celltype}_opt_fpc.csv",
        vcf="../../data/sample_info/fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz",
        bed="{celltype}_ca_75_normalized.bed.gz",
    output:
        "../../results/results_06032025/QTL_opt_results/{celltype}_{region}/nominal_chr1.txt"
    params:
        celltype="{celltype}",
        region="ca",
        nocov_dir="../../results/results_08312024/QTLscan_nocov_permutation_75ext/{celltype}",
        cov="../../results/results_08312024/covariates_75ext/qtl_cov_simple_{celltype}_ca_",
        outdir="/scratch/scjp_root/scjp1/xiaoouw/Muscle_snATAC_nucQTL/results/results_06032025/QTL_opt_results/{celltype}_ca/",
        window_size='500000',
        localcov="{celltype}_ca_simple_cov.txt",
        heatmap_dir="../../results/results_08312024/heatmaps_covariates_75ext",
        local_bed="{celltype}_ca_75_normalized.bed",
        vcf="fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz"
    envmodules:
        "singularity",
    shell:
        """
        scp {input.vcf}* .
        fpc=$(python -c "import pandas as pd; df = pd.read_csv('{input.summary}', index_col=0); print(df.loc['{params.celltype}', '{params.region}'])")
        echo ${{fpc}}
        scp {params.nocov_dir}_{params.region}_F${{fpc}}_G5/merged_QTLresults_fdrcorr.csv {params.outdir}merged_QTLresults_fdrcorr.csv
        # heatmaps
        scp {params.heatmap_dir}/FPC${{fpc}}_{params.celltype}_{params.region}.png {params.outdir}FPC${{fpc}}_{params.celltype}_{params.region}.png
        # get cov
        scp {params.cov}FPC${{fpc}}_GPC5.txt {params.outdir}{params.celltype}_{params.region}_simple_cov.txt
        # local cov
        scp {params.cov}FPC${{fpc}}_GPC5.txt {params.localcov}

        for i in 1
        do
            singularity exec /gpfs/accounts/scjp_root/scjp99/arushiv/muscle-sn/data/containers/qtl/qtl.sif QTLtools cis \
                    --vcf {params.vcf} \
                    --bed {params.local_bed}.gz \
                    --cov {params.localcov} \
                    --region chr${{i}} \
                    --nominal 1.0 \
                    --out {params.outdir}nominal_chr${{i}}.txt \
                    --normal \
                    --window {params.window_size} \
                    --std-err \
                    --seed 123
        done
        rm -r {params.localcov}*
        rm -r {params.local_bed}*
        """
rule correct_bed:
   input: 
       bed="{celltype}_ca_75_normalized.bed.gz",
       cov="../../results/results_06032025/QTL_opt_results/{celltype}_ca/{celltype}_ca_simple_cov.txt"
   params:
       localcov="{celltype}_ca_simple_cov.txt"
   envmodules:
         "singularity",
   output:
       "/scratch/scjp_root/scjp1/xiaoouw/Muscle_snATAC_nucQTL/results/results_06032025/corrected_bed/{celltype}_ca_corrected.bed.gz"
   shell:
       """
       scp {input.cov} {params.localcov}
       singularity exec /gpfs/accounts/scjp_root/scjp99/arushiv/muscle-sn/data/containers/qtl/qtl.sif QTLtools correct --bed {input.bed} \
                        --cov {params.localcov}\
                        --normal \
                        --out {output}
        rm -r {params.localcov}
       """

rule runsusie:
    input:
        bed="../../results/results_06032025/corrected_bed/{celltype}_ca_corrected.bed.gz",
        qtl="../../results/results_06032025/QTL_opt_results/{celltype}_ca/merged_QTLresults_fdrcorr.csv"       
    params:
        chr="chr{id}",
        min_corr=config["min_corr"],
        num_L=config["num_L"],
        outdir=config["outdir"],
        nominal_pass_dir=config["nominal_dir"],
        celltype="{celltype}",
        region="ca",
    output:
        "../../results/results_06032025/{celltype}_ca_susie/chr{id}/chr{id}_susie_cset95_summary.tsv"
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
            --nominal {params.nominal_pass_dir}{params.celltype}_{params.region}/nominal_{params.chr}.txt
        """


