configfile: "coloc_gwas_arushi.yaml"

traits=config["traits"]
filter_window=config["filter_window"]
PPH4_threshold=config["PPH4_threshold"]
celltypes=config["celltypes"]
regions=config["regions"]

import glob

# Function to list all files in a directory
def list_files_in_dir(directory, pattern="*"):
    return glob.glob(f"{directory}/{pattern}")

# Example: Listing all .txt files in the 'data' directory
input_files = list_files_in_dir("../../data/susie_gwas", "*.Rda")

rule all:
    input:
        expand("../../results/results_08152024/coloc_gwas_arushi/{file}_out.Rda", file=[f.replace("../../data/susie_gwas/", "").replace(".Rda", "_out.Rda") for f in input_files])

rule coloc:
    input:
        "../../data/susie_gwas/{file}.Rda"
    params:
        celltype="{celltype}",
        indir_qtl = config["indir_qtl"],
        indir_gwas = config["indir_gwas"],
        outdir_coloc = config["outdir_coloc"],
        filter_window=config["filter_window"],
        region="{region}",
    output:
        "../../results/results_08152024/coloc_gwas_arushi/{file}_out.Rda"
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        Rscript ../../scripts/bin/coloc_gwas_arushi.R \
            --gwas_dat {input} \
            --indir_qtl {params.indir_qtl}{params.celltype}_{params.region}_susie/ \
            --outdir_coloc {params.outdir_coloc} \
            --filter_window {params.filter_window} \
            --region {params.region}
        """