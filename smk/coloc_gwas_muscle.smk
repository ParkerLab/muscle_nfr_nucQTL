configfile: "coloc_gwas_muscle.yaml"
celltypes=config['celltypes']
traits=config['traits']
filter_window=config["filter_window"]

import glob


rule all:
    input:
        #expand("../../results/results_08222024/ukbb_dat/{trait}/{trait}_leadvar.tsv", trait=traits),
        #expand("../../results/results_08262024/{trait}_susie_models/{trait}_summary_susie_cset95.tsv", trait=traits),
        expand("../../results/results_08272024/GWAS_polymyositis_coloc/{celltype}_nfr/gwas_chr11-67929137-68929137_res_summary_{filter_window}win.Rda", 
                celltype=celltypes, filter_window=filter_window),
        expand("../../results/results_08272024/gwas_coloc_plots_{celltype}_{filter_window}win/polymyositis_nfr/locuscomapreplot/logs/gwas_chr11-67929137-68929137_win_PPH40.5.tsv",
                celltype=celltypes, filter_window=filter_window),

rule get_lead_var:
    input:
        '../../data/GWAS/{trait}.h.tsv.gz'
    output:
        "../../results/results_08222024/ukbb_dat/{trait}/{trait}_leadvar.tsv"
    params:
        trait="{trait}",
        outdir="../../results/results_08222024/ukbb_dat/{trait}/"
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        python ../../scripts/bin/get_gwas_lead_gwascat.py --summary {input} \
            --window 250000 \
            --traitname {params.trait} \
            --output {params.outdir}{params.trait}_leadvar.tsv
        """
rule runsusie:
    input:
        "../../results/results_08222024/ukbb_dat/{trait}/{trait}_leadvar.tsv"
    output:
        "../../results/results_08262024/{trait}_susie_models/{trait}_summary_susie_cset95.tsv"
    params:
        trait="{trait}",
        sample_size=lambda wildcards: config["sample_size"][wildcards.trait],
        ukbb_path="../../data/ukbb_hg38/",
        gwas_file_path="../../data/GWAS/{trait}.h.tsv.gz",
        outdir="../../results/results_08262024/"
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        Rscript ../../scripts/bin/runsusie_gwas_muscle.R --GWAS_annot_leadvar {input} \
            --gwas_file_path {params.gwas_file_path} \
            --sample_size {params.sample_size} \
            --ukbb_path {params.ukbb_path} \
            --outdir {params.outdir}
        """

rule coloc:
    input:
        "../../results/results_08262024/polymyositis_susie_models/gwas_chr11-67929137-68929137.susie.Rda",
    params:
        LD_file = "gwas_chr11-67929137-68929137",
        celltype="{celltype}",
        GWAS_annot="polymyositis",
        indir_qtl = "../../results/results_05132024/",
        indir_gwas = "../../results/results_08262024/",
        outdir_coloc = "../../results/results_08272024/",
        filter_window=config["filter_window"],
        region="nfr",
    output:
        "../../results/results_08272024/GWAS_polymyositis_coloc/{celltype}_nfr/gwas_chr11-67929137-68929137_res_summary_{filter_window}win.Rda"
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        Rscript ../../scripts/bin/coloc_gwas.R \
            --gwas_dat {input} \
            --GWAS_annot {params.GWAS_annot} \
            --indir_qtl {params.indir_qtl}{params.celltype}_{params.region}_susie/ \
            --indir_gwas {params.indir_gwas} \
            --outdir_coloc {params.outdir_coloc}GWAS_{params.GWAS_annot}_coloc/{params.celltype}_{params.region}/ \
            --filter_window {params.filter_window} \
            --region {params.region}
        """

rule tidy_coloc:
    input:
        "../../results/results_08272024/GWAS_polymyositis_coloc/{celltype}_nfr/gwas_chr11-67929137-68929137_res_summary_{filter_window}win.Rda"
    output:
        "../../results/results_08272024/gwas_coloc_plots_{celltype}_{filter_window}win/polymyositis_nfr/locuscomapreplot/logs/gwas_chr11-67929137-68929137_win_PPH40.5.tsv"
    params:
        LD_file = "gwas_chr11-67929137-68929137",
        celltype="{celltype}",
        outdir_coloc="../../results/results_08272024/",
        GWAS_annot="polymyositis",
        region="nfr",
        PPH4_threshold="0.5",
        dir_qtl_susie="../../results/results_05132024/",
        dir_gwas_susie="../../results/results_08262024/",
        filter_window=config["filter_window"],
        outdir_plot="../../results/results_08272024/",
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        Rscript ../../scripts/bin/tidy_gwas_coloc.R \
            --LD_file {params.LD_file} \
            --celltype {params.celltype} \
            --outdir_coloc {params.outdir_coloc} \
            --GWAS_annot {params.GWAS_annot} \
            --region {params.region} \
            --PPH4_threshold {params.PPH4_threshold} \
            --outdir_plot {params.outdir_plot} \
            --dir_qtl_susie {params.dir_qtl_susie} \
            --dir_gwas_susie {params.dir_gwas_susie} \
            --filter_window {params.filter_window} \
            --gwas_summary {input}
        #cat {input} > {output}
        """