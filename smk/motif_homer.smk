configfile: "motif_homer.yaml"

celltypes=config["celltypes"]

rule all:
    input:
        #expand("../../results/results_10032024/motif_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_nuc.homerpeak", celltype=celltypes),
        #expand("../../results/results_10032024/motif_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_nuc_bkg.homerpeak", celltype=celltypes),
        expand("../../results/results_02182025/motif_coloc/{celltype}_nfrqtl_coloc/homerResults.html", celltype=celltypes),
        #expand("../../results/results_10032024/motif_gwas_coloc/{celltype}_nfrqtl_coloc/t2d_homerResults.html", celltype=celltypes)

rule tify_coloc_motif:
    input:
        nfrqtl="../../results/results_05102024/QTL_opt_results/{celltype}_nfr/merged_QTLresults_fdrcorr.csv",
        coloc="../../results/results_07172024/coloc_plots_{celltype}/coloc_summary_{celltype}_5000win_PPH40.5.tsv",
    params:
        celltype="{celltype}",
        outdir="../../results/results_02182025/motif_coloc/{celltype}_nfrqtl_coloc/",
    output:
        "../../results/results_02182025/motif_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_nuc.homerpeak",
        "../../results/results_02182025/motif_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_nuc_bkg.homerpeak",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/tidy_nfr_coloc_motif.R \
            --nfrqtl {input.nfrqtl} \
            --coloc {input.coloc} \
            --celltype {params.celltype} \
            --outdir {params.outdir}
        """

rule coloc_motif:
    input:
        nfrpeak="../../results/results_02182025/motif_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_nuc.homerpeak",
        bgpeak="../../results/results_02182025/motif_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_nuc_bkg.homerpeak"
    params:
        celltype="{celltype}",
        outdir="../../results/results_02182025/motif_coloc/{celltype}_nfrqtl_coloc/",
    output:
        "../../results/results_02182025/motif_coloc/{celltype}_nfrqtl_coloc/homerResults.html"
    envmodules:
        "Bioinformatics",
        "homer/4.11.1-3coenyl",
    shell:
        """
        findMotifsGenome.pl \
            {input.nfrpeak} \
            hg38 \
            {params.outdir} \
            -bg {input.bgpeak} \
            -size given
        """

rule tify_gwas_motif:
    input:
        nfrqtl="../../results/results_05102024/QTL_opt_results/{celltype}_nfr/merged_QTLresults_fdrcorr.csv",
        gwas="../../results/results_05212024/nfrqtls_gwas_coloc_05142024.tsv",
    params:
        celltype="{celltype}",
        outdir="../../results/results_10032024/motif_gwas_coloc/{celltype}_nfrqtl_coloc/",
        trait="t2d",
    output:
        "../../results/results_10032024/motif_gwas_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_gwas_t2d.homerpeak",
        "../../results/results_10032024/motif_gwas_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_gwas_t2d_bkg.homerpeak",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/tidy_nfr_gwas_motif.R \
            --nfrqtl {input.nfrqtl} \
            --gwas {input.gwas} \
            --celltype {params.celltype} \
            --outdir {params.outdir} \
            --trait {params.trait}
        """

rule gwas_motif:
    input:
        nfrpeak="../../results/results_10032024/motif_gwas_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_gwas_t2d.homerpeak",
        bgpeak="../../results/results_10032024/motif_gwas_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_gwas_t2d_bkg.homerpeak"
    params:
        celltype="{celltype}",
        outdir="../../results/results_10032024/motif_gwas_coloc/{celltype}_nfrqtl_coloc/",
        trait="t2d",
    output:
        "../../results/results_10032024/motif_gwas_coloc/{celltype}_nfrqtl_coloc/t2d_homerResults.html"
    envmodules:
        "Bioinformatics",
        "homer/4.11.1-3coenyl",
    shell:
        """
        findMotifsGenome.pl \
            {input.nfrpeak} \
            hg38 \
            {params.outdir}{params.trait}_ \
            -bg {input.bgpeak} \
            -size given
        """