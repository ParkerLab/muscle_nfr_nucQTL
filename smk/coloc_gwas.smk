configfile: "coloc_gwas.yaml"
GWAS_annots=config["GWAS_annotations"]
ids=config["IDS"]
traits=config["traits"]
filter_window=config["filter_window"]
types=config["types"]
PPH4_threshold=config["PPH4_threshold"]
celltypes=config["celltypes"]
regions=config["regions"]

import pandas as pd

def get_LD_file(lead_var_file):
    lead_var_df = pd.read_csv(lead_var_file, sep = '\t')
    # init lists
    LD_files = []
    for i in range(len(lead_var_df)):
        LD_chr = str(lead_var_df.iloc[i,]['Chr'])
        LD_start = str(lead_var_df.iloc[i,]['Start'])
        LD_end = str(lead_var_df.iloc[i,]['End'])
        LD_file = 'gwas_chr' + LD_chr+ "-"+ LD_start+ "-"+ LD_end
        LD_files.append(LD_file)
    return(LD_files)

LD_files = get_LD_file("../../data/GWAS/diamante_T2D-European_gwas_lead_variants_hg38_ext_regions.txt")

rule all:
    input:
        #expand("../../results/results_06262023/{annotation}_sigQTL_susie_models/{annotation}_chr{id}_done.txt", annotation=config["annotations"].keys(), id=ids)
        #expand("../../results/results_07052023/coloc_{trait}/chr{id}/chr{id}_res_summary_{filter_window}win.Rda", trait=traits, id=ids, filter_window=filter_window)
        #expand("../../results/results_07132023/{GWAS_annot}_susie_models/chr{id}/chr{id}_susie_cset95_summary.tsv", GWAS_annot=GWAS_annots, id=ids)
        #expand("../../results/results_07142023/{GWAS_annot}_susie_models/{LD_file}.susie_pipdf.tsv", GWAS_annot=GWAS_annots, LD_file=LD_files)
        #expand("../../results/results_07142023/{GWAS_annot}_susie_models/{LD_file}_alignment.png", GWAS_annot=GWAS_annots, LD_file=LD_files),
        ##expand("../../results/results_05212024/GWAS_{GWAS_annot}_coloc/{celltype}_{region}/{LD_file}_res_summary_{filter_window}win.Rda", celltype=celltypes, region=regions, GWAS_annot=GWAS_annots, LD_file=LD_files, filter_window=filter_window)
        expand("../../results/results_05222024/gwas_coloc_plots_{celltype}_{filter_window}win/{GWAS_annot}_{region}/locuscomapreplot/logs/{LD_file}_win_PPH4{PPH4_threshold}.tsv", 
                celltype=celltypes,region=regions, GWAS_annot=GWAS_annots, LD_file=LD_files, PPH4_threshold=PPH4_threshold, filter_window=filter_window),
        expand('../../results/results_03312025/GWAS_{GWAS_annot}_coloc_summary/nfrqtls_{GWAS_annot}_gwas_coloc_{celltype}_chr{id}_inpeak.tsv', 
                GWAS_annot=GWAS_annots, celltype=celltypes, id=ids),
        expand('../../results/results_03312025/GWAS_{GWAS_annot}_coloc_inpeak/nfrqtls_{GWAS_annot}_gwas_coloc_{celltype}_inpeak.tsv', 
                GWAS_annot=GWAS_annots, celltype=celltypes)


rule runsusie:
    input:
        "../../data/GWAS/{GWAS_annot}.bed.gz"
    params:
        GWAS_annot="{GWAS_annot}",
        LD_file = "{LD_file}",
        ukbb_path = config["VCF_path"],
        num_L = config["num_L"],
        min_corr = config["min_corr"],
        max_it = config["max_it"],
        outdir = config["outdir_susie"]
    output:
        "../../results/results_07142023/{GWAS_annot}_susie_models/{LD_file}_susie_dat.Rda"
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        Rscript ../../scripts/bin/runsusie_gwas.R \
            --GWAS_annot {params.GWAS_annot} \
            --gwas_file_path {input} \
            --LD_file {params.LD_file} \
            --ukbb_path {params.ukbb_path} \
            --min_corr {params.min_corr} \
            --num_L {params.num_L} \
            --max_it {params.max_it} \
            --outdir {params.outdir}
        """

rule plot:
    input:
        "../../results/results_07142023/{GWAS_annot}_susie_models/{LD_file}_susie_dat.Rda"
    params:
        outdir_susie=config["outdir_susie"],
        GWAS_annot="{GWAS_annot}",
        LD_file = "{LD_file}",
    output:
        "../../results/results_07142023/{GWAS_annot}_susie_models/{LD_file}_alignment.png",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/plot_alignment.R \
            --filepath {input} \
            --out_filepath {params.outdir_susie}/{params.GWAS_annot}_susie_models/{params.LD_file}_alignment
        """

rule coloc:
    input:
        "../../results/results_07142023/{GWAS_annot}_susie_models/{LD_file}.susie.Rda",
    params:
        LD_file = "{LD_file}",
        celltype="{celltype}",
        GWAS_annot="{GWAS_annot}",
        indir_qtl = config["indir_qtl"],
        indir_gwas = config["indir_gwas"],
        outdir_coloc = config["outdir_coloc"],
        filter_window=config["filter_window"],
        region="{region}",
    output:
        "../../results/results_05212024/GWAS_{GWAS_annot}_coloc/{celltype}_{region}/{LD_file}_res_summary_{filter_window}win.Rda"
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
        "../../results/results_05212024/GWAS_{GWAS_annot}_coloc/{celltype}_{region}/{LD_file}_res_summary_{filter_window}win.Rda"
    output:
        "../../results/results_05222024/gwas_coloc_plots_{celltype}_{filter_window}win/{GWAS_annot}_{region}/locuscomapreplot/logs/{LD_file}_win_PPH4{PPH4_threshold}.tsv"
    params:
        LD_file = "{LD_file}",
        celltype="{celltype}",
        outdir_coloc=config["outdir_coloc"],
        GWAS_annot="{GWAS_annot}",
        region="{region}",
        PPH4_threshold=config["PPH4_threshold"],
        dir_qtl_susie=config["indir_qtl"],
        dir_gwas_susie=config["indir_gwas"],
        filter_window=config["filter_window"],
        outdir_plot=config["outdir_plot"],
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

rule coloc_summary:
    input:
        "../../data/GWAS/{GWAS_annot}.bed.gz"
    output:
        "../../results/results_03312025/GWAS_{GWAS_annot}_coloc_summary/nfrqtls_{GWAS_annot}_gwas_coloc.tsv"
    params:
        GWAS_annot="{GWAS_annot}",
        rda_dir = '../../results/results_05212024/GWAS_{GWAS_annot}_coloc/',
        filter_window=config["filter_window"],
        outdir= '../../results/results_03312025/GWAS_{GWAS_annot}_coloc_summary/'
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/coloc_summary.R \
            --rda_dir {params.rda_dir} \
            --filter_window {params.filter_window} \
            --GWAS_annot {params.GWAS_annot} \
            --outdir {params.outdir}
        """

rule coloc_inpeak:
    input:
        "../../results/results_03312025/GWAS_{GWAS_annot}_coloc_summary/nfrqtls_{GWAS_annot}_gwas_coloc.tsv"
    output:
        "../../results/results_03312025/GWAS_{GWAS_annot}_coloc_summary/nfrqtls_{GWAS_annot}_gwas_coloc_{celltype}_chr{id}_inpeak.tsv"
    params:
        GWAS_annot="{GWAS_annot}",
        nominal_path=config["nominal_path"],
        nfrqtl_susie_dir=config["indir_qtl"],
        gwas_susie_dir=config["indir_gwas"],
        chr="chr{id}",
        celltype="{celltype}"
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/coloc_inpeak_var.R \
            --celltype {params.celltype} \
            --chr {params.chr} \
            --nominal_path {params.nominal_path} \
            --nfrqtl_susie_dir {params.nfrqtl_susie_dir} \
            --gwas_susie_dir {params.gwas_susie_dir}/{params.GWAS_annot}_susie_models/ \
            --gwas_summary {input} \
            --output {output}
        """

rule coloc_inpeak_bind:
    input:
        "../../results/results_03312025/GWAS_{GWAS_annot}_coloc_summary/nfrqtls_{GWAS_annot}_gwas_coloc_{celltype}_chr22_inpeak.tsv"
    output:
        "../../results/results_03312025/GWAS_{GWAS_annot}_coloc_inpeak/nfrqtls_{GWAS_annot}_gwas_coloc_{celltype}_inpeak.tsv"
    params:
        GWAS_annot="{GWAS_annot}",
        celltype="{celltype}",
        inpeak_dir="../../results/results_03312025/GWAS_{GWAS_annot}_coloc_summary/",
        outdir="../../results/results_03312025/GWAS_{GWAS_annot}_coloc_inpeak/"
    shell:
        """
        # Define the pattern and output
        input_directory="{params.inpeak_dir}"
        pattern="*{params.celltype}*.tsv"
        output_file="{output}"

        # Initialize the output file
        > "$output_file"  # Create an empty output file

        # Check to include the header from the first non-empty file
        first=true

        for file in "$input_directory"/$pattern; do
        # Skip files that do not exist or are empty
        if [ -s "$file" ]; then  # -s checks if the file is not empty
            if $first; then
            # Add header from the first non-empty file
            head -n 1 "$file" >> "$output_file"
            first=false
            fi
            # Append the data from each file, skipping the header row
            tail -n +2 "$file" >> "$output_file"
        fi
        done
        """