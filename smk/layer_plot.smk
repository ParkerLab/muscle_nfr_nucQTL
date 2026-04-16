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
        motif_db = config.get("motif_db", "../../data/motif/H12CORE_pfms.txt"),
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
        nominal_path=config.get("nominal_path", "../../results/results_05102024/QTL_opt_results/"),
        vcf_file=config.get("vcf_file", "../../data/sample_info/fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz"),
        bigwig_sample_map_nfr=config.get("bigwig_sample_map_nfr", "../../data/bw_maps/{cluster}_nfr.map"),
        bigwig_sample_map_nuc=config.get("bigwig_sample_map_nuc", "../../data/bw_maps/{cluster}_nuc.map"),
        bigwig_sample_map_nuc_smooth=config.get("bigwig_sample_map_nuc_smooth", "../../data/bw_maps/{cluster}_nuc_smooth.map"),
        font_path=config.get("font_path", ""),
    conda:
        "pygenometracks"
    shell:
        """
        python ../../scripts/bin/plotlayered.py {input} {params.outdir} {params.nominal_path} {params.vcf_file} {params.bigwig_sample_map_nfr} {params.bigwig_sample_map_nuc} {params.bigwig_sample_map_nuc_smooth} {params.font_path}
        """
