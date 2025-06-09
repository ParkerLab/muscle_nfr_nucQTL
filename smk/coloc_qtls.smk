configfile: "coloc_qtls.yaml"

ids=config["IDS"]
annotations=list(config["annotations"].keys())
traits=config["traits"]
filter_window=config["filter_window"]
PPH4_threshold=config["PPH4_threshold"]

celltypes=config["celltypes"]
regions=config["regions"]

# wildcard_constraints:
#     annotation="|".join(annotations),

# def get_method(annotation):
#     return(config["annotations"][annotation])


rule all:
    input:
        #expand("../../data/QTL_phenotype_data/temp_beds/{annotation}_corrected.Rda", annotation=config["annotations"].keys())
        #expand("../../results/results_05122024/corrected_bed/{celltype}_{region}_corrected.bed.gz", celltype=celltypes, region=regions),
        #expand("../../results/results_07162024/{celltype}_{region}_susie/chr{id}/chr{id}_susie_cset95_summary.tsv", celltype=celltypes, region=regions, id=ids),
        #expand("../../results/results_07172024/coloc_{celltype}/chr{id}/chr{id}_res_summary_{filter_window}win.Rda", celltype=celltypes, id=ids, filter_window=filter_window),
        #expand("../../results/results_10172023/{annotation}_sigQTL_susie_models/chr{id}/chr{id}_susie_cset95_summary.tsv", annotation=config["annotations"].keys(), id=ids)
        #expand("../../results/results_06262023/{annotation}_sigQTL_susie_models/{annotation}_chr{id}_done.txt", annotation=config["annotations"].keys(), id=ids)
        #expand("../../results/results_09112023/coloc_{trait}/chr{id}/chr{id}_res_summary_{filter_window}win.Rda", trait=traits, id=ids, filter_window=filter_window)
        expand("../../results/results_07172024/coloc_plots_{celltype}/coloc_summary_{celltype}_{filter_window}win_PPH4{PPH4_threshold}.tsv", celltype=celltypes, filter_window=filter_window, PPH4_threshold=PPH4_threshold),
        expand("../../results/results_09102024/coloc_regression_plot/coloc_regression_plot_{celltype}.tsv", celltype=celltypes)

rule correct_bed:
   input: 
       bed="{celltype}_{region}_normalized.bed.gz",
       cov="../../results/results_07162024/QTL_opt_results/{celltype}_{region}/{celltype}_{region}_simple_cov.txt"
   params:
       localcov="{celltype}_{region}_simple_cov.txt"
   output:
       "../../results/results_07162024/corrected_bed/{celltype}_{region}_corrected.bed.gz"
   shell:
       """
       scp {input.cov} {params.localcov}
       QTLtools correct --bed {input.bed} \
                        --cov {params.localcov}\
                        --normal \
                        --out {output}
        rm -r {params.localcov}
       """

rule runsusie:
    input:
        bed="../../results/results_07162024/corrected_bed/{celltype}_{region}_corrected.bed.gz",
        qtl="../../results/results_07162024/QTL_opt_results/{celltype}_{region}/merged_QTLresults_fdrcorr.csv"       
    params:
        chr="chr{id}",
        min_corr=config["min_corr"],
        num_L=config["num_L"],
        outdir=config["outdir"],
        nominal_pass_dir=config["nominal_dir"],
        celltype="{celltype}",
        region="{region}",
    output:
        "../../results/results_07162024/{celltype}_{region}_susie/chr{id}/chr{id}_susie_cset95_summary.tsv"
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

rule coloc:
    input:
        nfr="../../results/results_05132024/{celltype}_nfr_susie/chr{id}/chr{id}_susie_cset95_summary.tsv",
        nuc="../../results/results_07162024/{celltype}_nuc_susie/chr{id}/chr{id}_susie_cset95_summary.tsv",
    output:
        "../../results/results_07172024/coloc_{celltype}/chr{id}/chr{id}_res_summary_{filter_window}win.Rda"
    params:
        chr="chr{id}",
        celltype="{celltype}",
        indir_nuc="../../results/results_07162024/",
        indir_nfr="../../results/results_05132024/",
        outdir_coloc=config["outdir_coloc"],
        filter_window=config["filter_window"],
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        Rscript ../../scripts/bin/coloc.R \
            --chr {params.chr} \
            --nfr {input.nfr} \
            --nuc {input.nuc} \
            --indir_nfr {params.indir_nfr}{params.celltype}_nfr_susie/{params.chr}/ \
            --indir_nuc {params.indir_nuc}{params.celltype}_nuc_susie/{params.chr}/ \
            --outdir_coloc {params.outdir_coloc}coloc_{params.celltype}/{params.chr}/ \
            --filter_window {params.filter_window}
        """

rule tidy_coloc:
    input:
        expand("../../results/results_07172024/coloc_{celltype}/chr{id}/chr{id}_res_summary_{filter_window}win.Rda", id=ids, celltype=celltypes, filter_window=filter_window)
    output:
        "../../results/results_07172024/coloc_plots_{celltype}/coloc_summary_{celltype}_{filter_window}win_PPH4{PPH4_threshold}.tsv"
    params:
        celltype="{celltype}",
        outdir_coloc=config["outdir_coloc"],
        outdir_plot=config["outdir_plot"],
        outdir_susie=config["outdir"],
        filter_window=config["filter_window"],
        PPH4_threshold=config["PPH4_threshold"],
        count_mat_dir=config["count_mat_dir"],
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        mkdir -p {params.outdir_plot}coloc_plots_{params.celltype}/locuscomapreplot
        mkdir -p {params.outdir_plot}coloc_plots_{params.celltype}/scatterplot
        
        Rscript ../../scripts/bin/tidy_qtl_coloc.R \
            --celltype {params.celltype} \
            --outdir_coloc {params.outdir_coloc} \
            --outdir_plot {params.outdir_plot} \
            --outdir_susie {params.outdir_susie} \
            --PPH4_threshold {params.PPH4_threshold} \
            --filter_window {params.filter_window} \
            --nfr_count {params.count_mat_dir}{params.celltype}_nfr_normalized.rds \
            --nuc_count {params.count_mat_dir}{params.celltype}_nuc_normalized.rds
        """

rule coloc_regression:
    input:
        "../../results/results_07172024/nfr_nucqlts_coloc_07172024.tsv"
    output:
        "../../results/results_09102024/coloc_regression_plot/coloc_regression_plot_{celltype}.png",
        "../../results/results_09102024/coloc_regression_plot/coloc_regression_plot_{celltype}.tsv",
    params:
        celltype="{celltype}",
        output_name="../../results/results_09102024/coloc_regression_plot/coloc_regression_plot_{celltype}"
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/coloc_regression_plot.R \
            --celltype {params.celltype} \
            --coloc_result {input} \
            --output {params.output_name}
        """


