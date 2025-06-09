configfile: "coloc_gwas_ukbb.yaml"
celltypes=config['celltypes']
traits=config['traits']
filter_window=config["filter_window"]

import glob

rule all:
    input:
        expand("../../results/results_09232024/gwas_ukbb/{trait}/{trait}_hg38.h.tsv", trait=traits),
        expand("../../results/results_09232024/gwas_ukbb/{trait}/{trait}_leadvar.tsv", trait=traits),
        #expand("../../results/results_10072024/{trait}_susie_models/{trait}_summary_susie_cset95.tsv", trait=traits),
        #expand("../../results/results_08272024/GWAS_polymyositis_coloc/{celltype}_nfr/gwas_chr11-67929137-68929137_res_summary_{filter_window}win.Rda", 
        #        celltype=celltypes, filter_window=filter_window),
        #expand("../../results/results_08272024/gwas_coloc_plots_{celltype}_{filter_window}win/polymyositis_nfr/locuscomapreplot/logs/gwas_chr11-67929137-68929137_win_PPH40.5.tsv",
        #        celltype=celltypes, filter_window=filter_window),

rule liftover:
    input:
        "../../data/GWAS_UKBB/{trait}.tsv.bgz"
    output:
        "../../results/results_09232024/gwas_ukbb/{trait}/{trait}_hg38.h.tsv"
    params:
        vcf="../../data/sample_info/fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/liftover.py --summary {input} \
            --vcf {params.vcf} \
            --output {output}
        """

rule get_lead_var:
    input:
        '../../results/results_09232024/gwas_ukbb/{trait}/{trait}_hg38.h.tsv'
    output:
        "../../results/results_09232024/gwas_ukbb/{trait}/{trait}_leadvar.tsv"
    params:
        trait="{trait}",
        outdir="../../results/results_08282024/gwas_ukbb/{trait}/"
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
            --output {output}
        """

rule runsusie:
    input:
        "../../results/results_09232024/gwas_ukbb/{trait}/{trait}_leadvar.tsv"
    output:
        "../../results/results_10072024/{trait}_susie_models/{trait}_summary_susie_cset95.tsv"
    params:
        trait="{trait}",
        sample_size=lambda wildcards: config["sample_size"][wildcards.trait],
        ukbb_path="../../data/ukbb_hg38/",
        gwas_file_path="../../results/results_09232024/gwas_ukbb/{trait}/{trait}_hg38.h.tsv",
        outdir="../../results/results_10072024/"
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
