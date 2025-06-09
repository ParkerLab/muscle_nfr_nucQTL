configfile: "layer_plot.yaml"

#vars=config['vars']
gwas_df=config['gwas_df']
inpeak_df=config['inpeak_df']
traits=config['gwas_traits']

rule all:
    input:
        expand(config['motif_dir']+"{trait}"+"_inpeak_motif/results_pval01.rds", trait=traits),
        expand(config["layer_plot_dir"]+"{trait}"+"/last.png", trait=traits),

rule motif:
    input:
        config['inpeak_df']+"{trait}"+"_snps_inpeak_info/query_snps_df.txt"
    output:
        config['motif_dir']+"{trait}"+"_inpeak_motif/results_pval01.rds"
    params:
        motif_db = "../../data/motif/H12CORE_pfms.txt",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/motif_mb.R \
            --inpeak_df {input} \
            --output {output} \
            --motif_db {params.motif_db}

        """
rule layer_plot:
    input:
        config['inpeak_df']+"{trait}"+"_snps_inpeak_info/list_layerplot_type2a.txt"
    output:
        config["layer_plot_dir"]+"{trait}"+"/last.png",
    params:
        outdir=config["layer_plot_dir"]+"{trait}"+"/",
    conda:
        "pygenometracks"
    shell:
        """
        python ../../scripts/bin/plotlayered.py {input} {params.outdir}
        """