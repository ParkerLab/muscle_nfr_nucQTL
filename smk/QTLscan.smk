configfile: "QTLscan.yaml"

celltypes=config['celltypes']
regions=config['regions']
fpc_nums=config["fpc_num"]
gpc_nums=config["gpc_num"]

rule all:
    input:
        #expand("../../results/results_07152024/covariates/qtl_cov_full_{celltype}_{region}_FPC{fpc}_GPC{gpc}.txt", celltype=celltypes, region=regions, fpc=fpc_nums, gpc=gpc_nums),
        #expand("{celltype}_{region}_normalized.bed.gz", celltype=celltypes, region=regions),
        #expand("../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/permutation_chr22.txt", celltype=celltypes, region=regions, fpc=fpc_nums, gpc=gpc_nums),
        #expand("../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults.txt", celltype=celltypes, region=regions, fpc=fpc_nums, gpc=gpc_nums),
        #expand("../../results/results_07152024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/permutation_chr22.txt", celltype=celltypes, region=regions, fpc=fpc_nums, gpc=gpc_nums),
        #expand("../../results/results_07152024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults.txt", celltype=celltypes, region=regions, fpc=fpc_nums, gpc=gpc_nums),
        #expand("../../results/results_05022024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults.txt", celltype=celltypes, region=regions, fpc=fpc_nums, gpc=gpc_nums),
        #expand("../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults_fdrcorr.csv", celltype=celltypes, region=regions, fpc=fpc_nums, gpc=gpc_nums),
        #expand("../../results/results_07152024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults_fdrcorr.csv", celltype=celltypes, region=regions, fpc=fpc_nums, gpc=gpc_nums),
        #expand("../../results/results_07152024/QTL_summary/{celltype}_qtl_by_fpcnum_summary.txt", celltype=celltypes),
        #expand("../../results/results_07152024/QTL_summary/{celltype}_opt_fpc.csv", celltype=celltypes)
        #expand("../../results/results_07162024/QTL_opt_results/{celltype}_{region}/nominal_chr22.txt", celltype=celltypes, region=regions)
        expand("../../results/results_06112024/QTL_opt_results/{celltype}_{region}_byfull/nominal_chr22.txt", celltype=celltypes, region=regions)
#rule create_outdir:
#    params:
#        outdir="/lab/work/xiaoouw/snATAC_m_QTLscan/data/results/{sample}_F{fpc}_G{gpc}/"
#    output:
#        "/lab/work/xiaoouw/snATAC_m_QTLscan/data/results/{sample}_F{fpc}_G{gpc}/here.txt"
#    shell:
#        """
#        mkdir -p {params.outdir} && touch {output}
#        """

rule get_covariates:
    input:
        #"../../results/results_05082024/normalized_counts/{celltype}_{region}_normalized.rds"
        #"../../results/results_06042024/normalized_counts_byfullbam/{celltype}_{region}_normalized.rds"
        "../../results/results_07152024/normalized_counts/{celltype}_{region}_normalized.rds"
    params:
        outdir="../../results/results_07152024/covariates/",
        fpc="{fpc}",
        gpc="{gpc}",
        vcf_eigenvec=config["vcf_eigenvec"],
        celltype="{celltype}",
        region = "{region}",
        figdir="../../results/results_07152024/heatmaps_covariates/",
        sample_info="../../data/sample_info/{celltype}_sample_level_covariates_atac.tsv"
    output:
        "../../results/results_07152024/covariates/qtl_cov_full_{celltype}_{region}_FPC{fpc}_GPC{gpc}.txt"
    conda:
        "atac"
    shell:
        """
        mkdir -p {params.figdir}
        python ../../scripts/bin/get_covariates.py {params.fpc} {params.gpc} {params.outdir} {params.figdir} {params.sample_info} {params.vcf_eigenvec} {input} {params.celltype} {params.region}
        """

rule preparebed:
    input:
        #bed="../../results/results_05082024/QTL_files/counts_bed/{celltype}_{region}_normalized.bed",
        #bed="../../results/results_06042024/QTL_files/counts_bed_byfullbam/{celltype}_{region}_normalized.bed",
        bed="../../results/results_07152024/QTL_files/counts_bed/{celltype}_{region}_normalized.bed"
    params:
        celltype="{celltype}",
        region = "{region}",
        localbed="{celltype}_{region}_normalized.bed",
    output:
        "{celltype}_{region}_normalized.bed.gz"
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
        bed="{celltype}_{region}_normalized.bed.gz",
        cov="../../results/results_07152024/covariates/qtl_cov_simple_{celltype}_{region}_FPC{fpc}_GPC{gpc}.txt",
    params:
        window_size=config["window_size"],
        fpc="{fpc}",
        gpc="{gpc}",
        celltype="{celltype}",
        region = "{region}",
        outdir="../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/",
        localcov="qtl_cov_simple_{celltype}_{region}_FPC{fpc}_GPC{gpc}.txt",
    output:
        "../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/permutation_chr22.txt"
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

rule QTL_scan_cov:
    input:
        vcf="fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz",
        bed="{celltype}_{region}_normalized.bed.gz",
        cov="../../results/results_07152024/covariates/qtl_cov_full_{celltype}_{region}_FPC{fpc}_GPC{gpc}.txt",
    params:
        window_size=config["window_size"],
        fpc="{fpc}",
        gpc="{gpc}",
        celltype="{celltype}",
        region = "{region}",
        outdir="../../results/results_07152024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/",
        localcov="qtl_cov_full_{celltype}_{region}_FPC{fpc}_GPC{gpc}.txt",
    output:
        "../../results/results_07152024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/permutation_chr22.txt"
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
        nocov="../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/permutation_chr22.txt",
        wcov="../../results/results_07152024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/permutation_chr22.txt"
    params:
        nocov_permutations="../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}",
        wcov_permutations="../../results/results_07152024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}"
    output:
        #"../../results/results_05082024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults.txt",
        #"../../results/results_05082024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults.txt"
        nocov="../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults.txt",
        wcov="../../results/results_07152024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults.txt"
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/Tidy_QTLresults.py {params.nocov_permutations} 0.05
        python ../../scripts/bin/Tidy_QTLresults.py {params.wcov_permutations} 0.05
        """


rule calc_qvalue:
    input:
        nocov="../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults.txt",
        wcov="../../results/results_07152024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults.txt"
        #"../../results/results_06112024/QTLscan_nocov_permutation_byfull/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults.txt"
    output:
        "../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults_fdrcorr.csv",
        "../../results/results_07152024/QTLscan_cov_permutation/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults_fdrcorr.csv"
        #"../../results/results_06112024/QTLscan_nocov_permutation_byfull/{celltype}_{region}_F{fpc}_G{gpc}/merged_QTLresults_fdrcorr.csv"
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/qvalue.R {input.nocov}
        Rscript ../../scripts/bin/qvalue.R {input.wcov}
        """

rule permutation_summary:
    input:
        nocov="../../results/results_07152024/QTLscan_nocov_permutation/{celltype}_nuc_F100_G5/merged_QTLresults_fdrcorr.csv",
        wcov="../../results/results_07152024/QTLscan_cov_permutation/{celltype}_nuc_F100_G5/merged_QTLresults_fdrcorr.csv"
    output:
        "../../results/results_07152024/QTL_summary/{celltype}_qtl_by_fpcnum_summary.txt"
    params:
        outdir="../../results/results_07152024/QTL_summary/",
        celltype="{celltype}",
        nocov_dir_nuc = "../../results/results_07152024/QTLscan_nocov_permutation/",
        wcov_dir_nuc = "../../results/results_07152024/QTLscan_cov_permutation/",
        nocov_dir_nfr = "../../results/results_05082024/QTLscan_nocov_permutation/",
        wcov_dir_nfr = "../../results/results_05082024/QTLscan_cov_permutation/",
    conda:
        "atac"
    shell:
        """
        # args: 1. celltype 2. fpc 3. nocov_dir 4. withcov_dir 5. outdir
        python ../bin/qtl_byfpc_summary.py {params.celltype} {params.nocov_dir_nfr} {params.wcov_dir_nfr} {params.nocov_dir_nuc} {params.wcov_dir_nuc} {params.outdir}
        """

rule select_fpc:
    input:
        "../../results/results_07152024/QTL_summary/{celltype}_qtl_by_fpcnum_summary.txt"
    output:
        "../../results/results_07152024/QTL_summary/{celltype}_opt_fpc.csv"
    params:
        celltype="{celltype}",
    conda:
        "atac"
    shell:
        """
        # args: 1. summary txt 2. celltype 3. output
        python ../bin/find_maxQTL.py {input} {params.celltype} {output}
        """

rule nominal_pass:
    input:
        summary="../../results/results_07152024/QTL_summary/{celltype}_opt_fpc.csv",
        vcf="fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz",
        bed="{celltype}_{region}_normalized.bed.gz",
    output:
        "../../results/results_07162024/QTL_opt_results/{celltype}_{region}/nominal_chr22.txt"
    params:
        celltype="{celltype}",
        region="{region}",
        nocov_dir="../../results/results_07152024/QTLscan_nocov_permutation/{celltype}",
        cov="../../results/results_07152024/covariates/qtl_cov_simple_{celltype}_{region}_",
        outdir="../../results/results_07162024/QTL_opt_results/{celltype}_{region}/",
        window_size='500000',
        localcov="{celltype}_{region}_simple_cov.txt",
        heatmap_dir="../../results/results_07152024/heatmaps_covariates",
    conda:
        "atac"
    shell:
        """
        fpc=$(python -c "import pandas as pd; df = pd.read_csv('{input.summary}', index_col=0); print(df.loc['{params.celltype}', '{params.region}'])")
        echo ${{fpc}}
        scp {params.nocov_dir}_{params.region}_F${{fpc}}_G5/merged_QTLresults_fdrcorr.csv {params.outdir}merged_QTLresults_fdrcorr.csv
        # heatmaps
        scp {params.heatmap_dir}/FPC${{fpc}}_{params.celltype}_{params.region}.png {params.outdir}FPC${{fpc}}_{params.celltype}_{params.region}.png
        # get cov
        scp {params.cov}FPC${{fpc}}_GPC5.txt {params.outdir}{params.celltype}_{params.region}_simple_cov.txt
        # local cov
        scp {params.cov}FPC${{fpc}}_GPC5.txt {params.localcov}

        for i in $(seq 1 22)
        do
            QTLtools cis \
                    --vcf {input.vcf} \
                    --bed {input.bed} \
                    --cov {params.localcov} \
                    --region chr${{i}} \
                    --nominal 1.0 \
                    --out {params.outdir}nominal_chr${{i}}.txt \
                    --normal \
                    --window {params.window_size} \
                    --std-err \
                    --seed 123
        done
        rm -r {params.localcov}
        """

rule nominal_pass_full:
    input:
        summary="../../results/results_07152024/QTL_summary/{celltype}_opt_fpc.csv",
        vcf="fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz",
        bed="{celltype}_{region}_full_normalized.bed.gz",
    output:
        "../../results/results_06112024/QTL_opt_results/{celltype}_{region}_byfull/nominal_chr22.txt"
    params:
        celltype="{celltype}",
        region="{region}",
        nocov_dir="../../results/results_06112024/QTLscan_nocov_permutation_byfull/{celltype}",
        cov="../../results/results_06112024/covariates/qtl_cov_simple_{celltype}_{region}_",
        outdir="../../results/results_06112024/QTL_opt_results/{celltype}_{region}_byfull/",
        window_size='500000',
        localcov="{celltype}_{region}_simple_cov_byfull.txt",
        #heatmap_dir="../../results/results_07152024/heatmaps_covariates",
    conda:
        "atac"
    shell:
        """
        fpc=$(python -c "import pandas as pd; df = pd.read_csv('{input.summary}', index_col=0); print(df.loc['{params.celltype}', '{params.region}'])")
        echo ${{fpc}}
        scp {params.nocov_dir}_{params.region}_F${{fpc}}_G5/merged_QTLresults_fdrcorr.csv {params.outdir}merged_QTLresults_fdrcorr.csv
        # get cov
        scp {params.cov}FPC${{fpc}}_GPC5.txt {params.outdir}{params.celltype}_{params.region}_simple_cov.txt
        # local cov
        scp {params.cov}FPC${{fpc}}_GPC5.txt {params.localcov}

        for i in $(seq 1 22)
        do
            QTLtools cis \
                    --vcf {input.vcf} \
                    --bed {input.bed} \
                    --cov {params.localcov} \
                    --region chr${{i}} \
                    --nominal 1.0 \
                    --out {params.outdir}nominal_chr${{i}}.txt \
                    --normal \
                    --window {params.window_size} \
                    --std-err \
                    --seed 123
        done
        rm -r {params.localcov}
        """