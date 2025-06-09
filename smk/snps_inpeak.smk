configfile: "snps_inpeak.yaml"

ids=config["IDS"]
celltypes=config["celltypes"]
traits=config["traits"]

rule all:
    input:
        expand("../../results/results_04232025/{celltype}_snps_inpeak_info/{celltype}_chr{id}_info.tsv", celltype=celltypes, id=ids),
        #expand("../../results/results_04082025/{trait}_snps_inpeak_info/{celltype}_chr{id}_info_gwas.tsv", celltype=celltypes, id=ids, trait=traits),
rule snps_inpeak:
    input:
        "../../results/results_07172024/nfr_nucqlts_coloc_07172024.tsv"
    output:
        "../../results/results_04232025/{celltype}_snps_inpeak_info/{celltype}_chr{id}_info.tsv"
    params:
        chr="chr{id}",
        celltype="{celltype}",
        susie_dir="../../results/results_07162024/",
        nfr_nominaldir="../../results/results_05102024/QTL_opt_results/",
        nuc_nominaldir="../../results/results_07162024/QTL_opt_results/"
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/var_inpeak.R \
            --coloc_dir {input} \
            --susiedir {params.susie_dir} \
            --chr {params.chr} \
            --celltype {params.celltype} \
            --nfr_nominaldir {params.nfr_nominaldir} \
            --nuc_nominaldir {params.nuc_nominaldir} \
            --outdir {output}

        """

rule snps_inpeak_gwas:
    input:
        summary=lambda wildcards: config["traits"][wildcards.trait]["summary"]
    output:
        "../../results/results_04082025/{trait}_snps_inpeak_info/{celltype}_chr{id}_info_gwas.tsv"
    params:
        chr="chr{id}",
        celltype="{celltype}",
        susie_dir="../../results/results_07162024/",
        nfr_nominaldir="../../results/results_05102024/QTL_opt_results/",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/var_inpeak_gwas.R \
            --gwas_summary {input.summary} \
            --susiedir {params.susie_dir} \
            --chr {params.chr} \
            --celltype {params.celltype} \
            --nfr_nominaldir {params.nfr_nominaldir} \
            --output {output}
        """""