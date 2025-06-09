configfile: "cit.yaml"

filter_window=config["filter_window"]
PPH4_threshold=config["PPH4_threshold"]
#traits=config["traits"]
celltypes=config['celltypes']


import pandas as pd

def get_pairs(trait, PPH4_threshold):
    trait_summary_path = "../../results/results_07172024/coloc_plots_"
    trait_summary_path += trait + "/coloc_summary_" + trait + "_5000.tsv"
    trait_summary_df = pd.read_csv(trait_summary_path, sep = '\t', header=0)
    # filter by PPH4
    trait_summary_df = trait_summary_df[trait_summary_df['PP.H4.abf'] > PPH4_threshold]
    # init lists
    cit_pairs = []
    for i in range(len(trait_summary_df)):
        hit1 = str(trait_summary_df.iloc[i,]['hit1'])
        hit2 = str(trait_summary_df.iloc[i,]['hit2'])
        NFR_peak = str(trait_summary_df.iloc[i,]['NFR_peak'])
        Nuc_peak = str(trait_summary_df.iloc[i,]['Nucleosomal_peak'])
        # paste all information in
        cit_pair = trait + '-' + hit1 + '-' + hit2 + '-'
        cit_pair += NFR_peak + '-' + Nuc_peak
        # append to lsit
        cit_pairs.append(cit_pair)
    return(cit_pairs)

# cit_pairs_HMMRATAC = get_pairs('HMMRATAC', PPH4_threshold)
# cit_pairs_NucleoATAC = get_pairs('NucleoATAC', PPH4_threshold)
# cit_pairs = []
# for item in [cit_pairs_HMMRATAC, cit_pairs_NucleoATAC]:
#     # appending elements to the flat_list
#     cit_pairs += item

# Get the pairs for each celltype
all_pairs = {celltype: get_pairs(celltype, PPH4_threshold) for celltype in celltypes}

# Create a list of all pair files for the rule input
cit_pairs = [pair for pairs in all_pairs.values() for pair in pairs]

#cit_pairs = cit_pairs

rule all:
    input:
        #expand("../../results/results_05222024/cit_PPH4_threshold{PPH4_threshold}/{cit_pair}_cit_perm_cov.Rda", PPH4_threshold=PPH4_threshold, cit_pair=cit_pairs)
        expand("../../results/results_07212024/{celltype}_citfdr_cov_{PPH4_threshold}_H4_summary.tsv", celltype=celltypes, PPH4_threshold=PPH4_threshold)

rule cit:
    input:
        expand("../../results/results_07172024/coloc_plots_{celltype}/coloc_summary_{celltype}_5000.tsv", celltype=celltypes)
    output:
        "../../results/results_07212024/cit_PPH4_threshold{PPH4_threshold}/{cit_pair}_cit_perm_cov.Rda"
    params:
        PPH4_threshold=config["PPH4_threshold"],
        outdir = config["outdir_cit"],
        #trait = lambda wildcards: '{trait}',
        #celltype = "{celltype}",
        pair = "{cit_pair}",
        cov = config["cov_dir"],
        count_mat_dir=config["count_mat_dir"],
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n"
    shell:
        """
        Rscript ../../scripts/bin/cit.R \
            --pair {params.pair} \
            --outdir {params.outdir} \
            --PPH4_threshold {params.PPH4_threshold} \
            --cov {params.cov} \
            --count_mat_dir {params.count_mat_dir}
        """

rule cit_fdr:
    input:
        expand("../../results/results_07212024/cit_PPH4_threshold{PPH4_threshold}/{cit_pair}_cit_perm_cov.Rda", cit_pair=cit_pairs, PPH4_threshold=PPH4_threshold)
    output:
        "../../results/results_07212024/{celltype}_citfdr_cov_{PPH4_threshold}_H4_summary.tsv"
    params:
        PPH4_threshold=config["PPH4_threshold"],
        cit_outdir=config["outdir_cit"],
        celltype="{celltype}",
        fdr_outdir=config["fdr_outdir"],
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "samtools/1.13-fwwss5n",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        Rscript ../../scripts/bin/cit_fdr.R \
            --celltype {params.celltype} \
            --cit_outdir {params.cit_outdir} \
            --PPH4_threshold {params.PPH4_threshold} \
            --fdr_outdir {params.fdr_outdir}
        """