configfile: "coloc_gwas_nonukbb.yaml"
celltypes=config['celltypes']
traits=config['traits']
filter_window=config["filter_window"]
regions=config["regions"]
ids=config["IDS"]

import pandas as pd
import glob
import os

def get_LD_file(lead_var_file):
    lead_var_df = pd.read_csv(lead_var_file, sep = '\t', names=['Chr', 'Start', 'End', 'Name'])
    # init lists
    LD_files = []
    # get the trait name
    trait = lead_var_file.split('/')[-2]
    for i in range(len(lead_var_df)):
        LD_chr = str(lead_var_df.iloc[i,]['Chr'])
        LD_start = str(lead_var_df.iloc[i,]['Start'])
        LD_end = str(lead_var_df.iloc[i,]['End'])
        LD_file = LD_chr+ "-"+ LD_start+ "-"+ LD_end + ":" + trait
        LD_files.append(LD_file)
    return(LD_files)

LD_files1 = get_LD_file("../../results/results_01162025/gwas_nonukbb/nielsen_atrial-fibrillation/nielsen_atrial-fibrillation_leadvar.tsv")
LD_files2 = get_LD_file("../../results/results_01162025/gwas_nonukbb/ukbb_creatinine/ukbb_creatinine_leadvar.tsv")
LD_files3 = get_LD_file("../../results/results_01162025/gwas_nonukbb/MAGIC_FI-EUR/MAGIC_FI-EUR_leadvar.tsv")
# merge lists
LD_files = LD_files1 + LD_files2 + LD_files3
LD_files = LD_files

def get_LD(files):
    file_names = [os.path.basename(f) for f in files]
    LDs = [filename.replace(".susie.Rda", "", 1) for filename in file_names]
    return(LDs)
naf_files = glob.glob(config["outdir_susie"]+"/nielsen_atrial-fibrillation_susie_models/*.susie.Rda")
creatinine_files = glob.glob(config["outdir_susie"]+"/ukbb_creatinine_susie_models/*.susie.Rda")
fi_files = glob.glob(config["outdir_susie"]+"/MAGIC_FI-EUR_susie_models/*.susie.Rda")

naf_LDs = get_LD(naf_files)
creatinine_LDs = get_LD(creatinine_files)
fi_LDs = get_LD(fi_files)
rule all:
    input:
        #expand("../../results/results_01162025/gwas_nonukbb/{trait}/{trait}_leadvar.tsv", trait=traits),
        #expand("../../results/results_01162025/track/{LD_file}.txt", LD_file=LD_files),
        #expand("../../results/results_03312025/GWAS_MAGIC_FI-EUR_coloc/{celltype}_nfr/{fi_LD}_res_summary_{filter_window}win.Rda", celltype=celltypes, fi_LD=fi_LDs, filter_window=filter_window),
        #expand("../../results/results_03312025/GWAS_nielsen_atrial-fibrillation_coloc/{celltype}_nfr/{naf_LD}_res_summary_{filter_window}win.Rda", celltype=celltypes, naf_LD=naf_LDs, filter_window=filter_window),
        #expand("../../results/results_03312025/GWAS_ukbb_creatinine_coloc/{celltype}_nfr/{creatinine_LD}_res_summary_{filter_window}win.Rda", celltype=celltypes, creatinine_LD=creatinine_LDs, filter_window=filter_window),
        expand("../../results/results_03312025/GWAS_{trait}_coloc_summary/nfrqtls_{trait}_gwas_coloc_{celltype}_chr{id}_inpeak.tsv", trait=traits, celltype=celltypes, id=ids),
        expand("../../results/results_03312025/GWAS_{trait}_coloc_inpeak/nfrqtls_{trait}_gwas_coloc_{celltype}_inpeak.tsv", trait=traits, celltype=celltypes),

rule get_lead_var:
    input:
        '../../data/GWAS/{trait}.bed.gz'
    output:
        "../../results/results_01162025/gwas_nonukbb/{trait}/{trait}_leadvar.tsv"
    params:
        trait="{trait}",
        outdir="../../results/results_01162025/gwas_nonukbb/{trait}/"
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        python ../../scripts/bin/get-leads.py --summary {input} \
            --window 250000 \
            --traitname {params.trait} \
            --output {output}
        """

rule runsusie:
    input:
        "../../results/results_01162025/gwas_nonukbb/MAGIC_FI-EUR/MAGIC_FI-EUR_leadvar.tsv"
    params:
        LD_file = "{LD_file}",
        ukbb_path = config["VCF_path"],
        num_L = config["num_L"],
        min_corr = config["min_corr"],
        max_it = config["max_it"],
        outdir = config["outdir_susie"]
    output:
        "../../results/results_01162025/track/{LD_file}.txt"
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        Rscript ../../scripts/bin/runsusie_nonukbb.R \
            --LD_file {params.LD_file} \
            --ukbb_path {params.ukbb_path} \
            --min_corr {params.min_corr} \
            --num_L {params.num_L} \
            --max_it {params.max_it} \
            --outdir {params.outdir}
        """

rule coloc_naf:
    input:
        '../../results/results_01162025/nielsen_atrial-fibrillation_susie_models/{naf_LD}.susie.Rda'
    params:
        naf_LD= "{naf_LD}",
        celltype="{celltype}",
        GWAS_annot="nielsen_atrial-fibrillation",
        indir_qtl = config["indir_qtl"],
        indir_gwas = config["indir_gwas"],
        outdir_coloc = config["outdir_coloc"],
        filter_window=config["filter_window"],
        region="nfr",
    output:
        "../../results/results_03312025/GWAS_nielsen_atrial-fibrillation_coloc/{celltype}_{region}/{naf_LD}_res_summary_{filter_window}win.Rda"
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
rule coloc_creatinine:
    input:
        '../../results/results_01162025/ukbb_creatinine_susie_models/{creatinine_LD}.susie.Rda'
    params:
        creatinine_LD= "{creatinine_LD}",
        celltype="{celltype}",
        GWAS_annot="ukbb_creatinine",
        indir_qtl = config["indir_qtl"],
        indir_gwas = config["indir_gwas"],
        outdir_coloc = config["outdir_coloc"],
        filter_window=config["filter_window"],
        region="nfr",
    output:
        "../../results/results_03312025/GWAS_ukbb_creatinine_coloc/{celltype}_{region}/{creatinine_LD}_res_summary_{filter_window}win.Rda"
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

rule coloc_fi:
    input:
        '../../results/results_01162025/MAGIC_FI-EUR_susie_models/{fi_LD}.susie.Rda'
    params:
        fi_LD= "{fi_LD}",
        celltype="{celltype}",
        GWAS_annot="MAGIC_FI-EUR",
        indir_qtl = config["indir_qtl"],
        indir_gwas = config["indir_gwas"],
        outdir_coloc = config["outdir_coloc"],
        filter_window=config["filter_window"],
        region="nfr",
    output:
        "../../results/results_03312025/GWAS_MAGIC_FI-EUR_coloc/{celltype}_nfr/{fi_LD}_res_summary_{filter_window}win.Rda"
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

rule coloc_summary:
    input:
        "../../data/GWAS/{trait}.bed.gz"
    output:
        "../../results/results_03312025/GWAS_{trait}_coloc_summary/nfrqtls_{trait}_gwas_coloc.tsv"
    params:
        GWAS_annot="{trait}",
        rda_dir = '../../results/results_03312025/GWAS_{trait}_coloc/',
        filter_window=config["filter_window"],
        outdir= '../../results/results_03312025/GWAS_{trait}_coloc_summary/'
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
        "../../results/results_03312025/GWAS_{trait}_coloc_summary/nfrqtls_{trait}_gwas_coloc.tsv"
    output:
        "../../results/results_03312025/GWAS_{trait}_coloc_summary/nfrqtls_{trait}_gwas_coloc_{celltype}_chr{id}_inpeak.tsv"
    params:
        GWAS_annot="{trait}",
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
        "../../results/results_03312025/GWAS_{trait}_coloc_summary/nfrqtls_{trait}_gwas_coloc_{celltype}_chr22_inpeak.tsv"
    output:
        "../../results/results_03312025/GWAS_{trait}_coloc_inpeak/nfrqtls_{trait}_gwas_coloc_{celltype}_inpeak.tsv"
    params:
        GWAS_annot="{trait}",
        celltype="{celltype}",
        inpeak_dir="../../results/results_03312025/GWAS_{trait}_coloc_summary/",
        outdir="../../results/results_03312025/GWAS_{trait}_coloc_inpeak/"
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