configfile: "eqtl_overlap.yaml"

celltypes=config["celltypes"]


rule all:
    input:
        expand("../../results/results_10162024/eqtl_overlap/{celltype}_eqtl_summary.tsv", celltype=celltypes),
        expand("../../results/results_10162024/eqtl_overlap/{celltype}_QTLoverlap.tsv", celltype=celltypes),
        expand("../../results/results_10162024/eqtl_overlap/{celltype}_QTL_fisher_exact.Rda", celltype=celltypes),
        expand("../../results/results_06052025/logit/{celltype}_QTL_logit.rds", celltype=celltypes)

rule eqtl_summary:
    input:
        "../../data/eqtl/fusion.Type_1--covlist.7.pheno-list.aa/fusion.Type_1--covlist.7.pheno-list.aa.ABCC4__L1-chr13_95314597_A_C-rs9302064.cset.bed"
    output:
        "../../results/results_10162024/eqtl_overlap/{celltype}_eqtl_summary.tsv"
    params:
        eqtl_dir="../../data/eqtl/",
        celltype="{celltype}",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/eqtl_summary.R --celltype {params.celltype} --eqtl_dir {params.eqtl_dir} --output {output}
        """

rule eqtl_overlap:
    input:
        "../../results/results_10162024/eqtl_overlap/{celltype}_eqtl_summary.tsv"
    output:
        "../../results/results_10162024/eqtl_overlap/{celltype}_QTLoverlap.tsv"
    params:
        celltype="{celltype}",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/eqtl_overlap.R --celltype {params.celltype} --eqtl {input} --output {output}
        """

rule fisher_exact:
    input:
        "../../results/results_10162024/eqtl_overlap/{celltype}_QTLoverlap.tsv"
    output:
        "../../results/results_10162024/eqtl_overlap/{celltype}_QTL_fisher_exact.Rda"
    params:
        celltype="{celltype}",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/fisher_exact.R --celltype {params.celltype} --overlap {input} --output {output}
        """

rule logit:
    input:
        "../../results/results_10162024/eqtl_overlap/{celltype}_QTLoverlap.tsv"
    output:
        "../../results/results_06052025/logit/{celltype}_QTL_logit.rds"
    params:
        celltype="{celltype}",
        nfr_susie_dir="../../results/results_07162024/",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/logit_eqtl_nuc.R \
            --celltype {params.celltype} \
            --input {input} \
            --nfr_susie_dir {params.nfr_susie_dir} \
            --outdir {output}
        """