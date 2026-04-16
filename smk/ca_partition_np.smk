configfile: "ca_partition_np.yaml"

celltypes = config["celltypes"]

OUTDIR = "../../results/results_03032026"
NFR_NUC_DIR = "../../results/results_10102024/motif_meme_qtl/peaks"


rule all:
    input:
        expand(OUTDIR + "/ca_partition_np_bp1/{celltype}/all_{celltype}.bed", celltype=celltypes),
        expand(OUTDIR + "/ca_partition_np_f50/{celltype}/all_{celltype}.bed", celltype=celltypes),
        expand(OUTDIR + "/upset_figs_bp1/{celltype}_upsetqtl.png", celltype=celltypes),
        expand(OUTDIR + "/upset_figs_f50/{celltype}_upsetqtl.png", celltype=celltypes),


rule format_caqtl_np:
    input:
        caqtl_csv="../../results/results_02262026/QTL_opt_results/{celltype}_ca/merged_QTLresults_fdrcorr.csv",
    output:
        caqtl_bed=OUTDIR + "/ca_partition_np_peaks/{celltype}_caqtl.bed",
    conda:
        "R423"
    shell:
        """
        mkdir -p {OUTDIR}/ca_partition_np_peaks
        Rscript ../../scripts/bin/qtl_csv_to_bed.R --csv {input.caqtl_csv} --out {output.caqtl_bed}
        """


rule intersect_bp1:
    input:
        nfr=NFR_NUC_DIR + "/{celltype}_nfrqtl.bed",
        nuc=NFR_NUC_DIR + "/{celltype}_nucqtl.bed",
        ca=OUTDIR + "/ca_partition_np_peaks/{celltype}_caqtl.bed",
    output:
        ca_nfr=OUTDIR + "/ca_partition_np_bp1/{celltype}/ca_nfr_{celltype}.bed",
        ca_nuc=OUTDIR + "/ca_partition_np_bp1/{celltype}/ca_nuc_{celltype}.bed",
        nfr_nuc=OUTDIR + "/ca_partition_np_bp1/{celltype}/nfr_nuc_{celltype}.bed",
        nfr_ca=OUTDIR + "/ca_partition_np_bp1/{celltype}/nfr_ca_{celltype}.bed",
        nuc_ca=OUTDIR + "/ca_partition_np_bp1/{celltype}/nuc_ca_{celltype}.bed",
        nuc_nfr=OUTDIR + "/ca_partition_np_bp1/{celltype}/nuc_nfr_{celltype}.bed",
        three_beds=OUTDIR + "/ca_partition_np_bp1/{celltype}/all_{celltype}.bed",
    params:
        temp=OUTDIR + "/ca_partition_np_bp1/{celltype}/nfr_nuc_precise_{celltype}.bed",
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm",
    shell:
        """
        mkdir -p {OUTDIR}/ca_partition_np_bp1/{wildcards.celltype}
        bedtools intersect -wa -a {input.ca} -b {input.nfr} > {output.ca_nfr}
        bedtools intersect -wa -a {input.ca} -b {input.nuc} > {output.ca_nuc}
        bedtools intersect -wa -a {input.nfr} -b {input.nuc} > {output.nfr_nuc}
        bedtools intersect -wa -a {input.nfr} -b {input.ca} > {output.nfr_ca}
        bedtools intersect -wa -a {input.nuc} -b {input.ca} > {output.nuc_ca}
        bedtools intersect -wa -a {input.nuc} -b {input.nfr} > {output.nuc_nfr}
        bedtools intersect -a {input.nfr} -b {input.nuc} > {params.temp}
        bedtools intersect -wa -a {input.ca} -b {params.temp} > {output.three_beds}
        """


rule intersect_f50:
    input:
        nfr=NFR_NUC_DIR + "/{celltype}_nfrqtl.bed",
        nuc=NFR_NUC_DIR + "/{celltype}_nucqtl.bed",
        ca=OUTDIR + "/ca_partition_np_peaks/{celltype}_caqtl.bed",
    output:
        ca_nfr=OUTDIR + "/ca_partition_np_f50/{celltype}/ca_nfr_{celltype}.bed",
        ca_nuc=OUTDIR + "/ca_partition_np_f50/{celltype}/ca_nuc_{celltype}.bed",
        nfr_nuc=OUTDIR + "/ca_partition_np_f50/{celltype}/nfr_nuc_{celltype}.bed",
        nfr_ca=OUTDIR + "/ca_partition_np_f50/{celltype}/nfr_ca_{celltype}.bed",
        nuc_ca=OUTDIR + "/ca_partition_np_f50/{celltype}/nuc_ca_{celltype}.bed",
        nuc_nfr=OUTDIR + "/ca_partition_np_f50/{celltype}/nuc_nfr_{celltype}.bed",
        three_beds=OUTDIR + "/ca_partition_np_f50/{celltype}/all_{celltype}.bed",
    params:
        temp=OUTDIR + "/ca_partition_np_f50/{celltype}/nfr_nuc_precise_{celltype}.bed",
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm",
    shell:
        """
        mkdir -p {OUTDIR}/ca_partition_np_f50/{wildcards.celltype}
        bedtools intersect -wa -a {input.ca} -b {input.nfr} -f 0.50 -r > {output.ca_nfr}
        bedtools intersect -wa -a {input.ca} -b {input.nuc} -f 0.50 -r > {output.ca_nuc}
        bedtools intersect -wa -a {input.nfr} -b {input.nuc} -f 0.50 -r > {output.nfr_nuc}
        bedtools intersect -wa -a {input.nfr} -b {input.ca} -f 0.50 -r > {output.nfr_ca}
        bedtools intersect -wa -a {input.nuc} -b {input.ca} -f 0.50 -r > {output.nuc_ca}
        bedtools intersect -wa -a {input.nuc} -b {input.nfr} -f 0.50 -r > {output.nuc_nfr}
        bedtools intersect -a {input.nfr} -b {input.nuc} -f 0.50 -r > {params.temp}
        bedtools intersect -wa -a {input.ca} -b {params.temp} -f 0.50 -r > {output.three_beds}
        """


rule upset_bp1:
    input:
        ca_nfr=OUTDIR + "/ca_partition_np_bp1/{celltype}/ca_nfr_{celltype}.bed",
        ca_nuc=OUTDIR + "/ca_partition_np_bp1/{celltype}/ca_nuc_{celltype}.bed",
        nfr_nuc=OUTDIR + "/ca_partition_np_bp1/{celltype}/nfr_nuc_{celltype}.bed",
        nfr_ca=OUTDIR + "/ca_partition_np_bp1/{celltype}/nfr_ca_{celltype}.bed",
        nuc_ca=OUTDIR + "/ca_partition_np_bp1/{celltype}/nuc_ca_{celltype}.bed",
        nuc_nfr=OUTDIR + "/ca_partition_np_bp1/{celltype}/nuc_nfr_{celltype}.bed",
        three_beds=OUTDIR + "/ca_partition_np_bp1/{celltype}/all_{celltype}.bed",
    output:
        OUTDIR + "/upset_figs_bp1/{celltype}_upsetqtl.png",
    params:
        celltype="{celltype}",
        nfr=NFR_NUC_DIR + "/{celltype}_nfrqtl.bed",
        nuc=NFR_NUC_DIR + "/{celltype}_nucqtl.bed",
        ca=OUTDIR + "/ca_partition_np_peaks/{celltype}_caqtl.bed",
    conda:
        "atac"
    shell:
        """
        mkdir -p {OUTDIR}/upset_figs_bp1
        python ../../scripts/bin/upset_qtl.py {input.ca_nfr} {input.ca_nuc} {input.nfr_nuc} {input.nfr_ca} \
                {input.nuc_ca} {input.nuc_nfr} {input.three_beds} \
                {params.celltype} {params.nfr} {params.nuc} {params.ca} {output}
        """


rule upset_f50:
    input:
        ca_nfr=OUTDIR + "/ca_partition_np_f50/{celltype}/ca_nfr_{celltype}.bed",
        ca_nuc=OUTDIR + "/ca_partition_np_f50/{celltype}/ca_nuc_{celltype}.bed",
        nfr_nuc=OUTDIR + "/ca_partition_np_f50/{celltype}/nfr_nuc_{celltype}.bed",
        nfr_ca=OUTDIR + "/ca_partition_np_f50/{celltype}/nfr_ca_{celltype}.bed",
        nuc_ca=OUTDIR + "/ca_partition_np_f50/{celltype}/nuc_ca_{celltype}.bed",
        nuc_nfr=OUTDIR + "/ca_partition_np_f50/{celltype}/nuc_nfr_{celltype}.bed",
        three_beds=OUTDIR + "/ca_partition_np_f50/{celltype}/all_{celltype}.bed",
    output:
        OUTDIR + "/upset_figs_f50/{celltype}_upsetqtl.png",
    params:
        celltype="{celltype}",
        nfr=NFR_NUC_DIR + "/{celltype}_nfrqtl.bed",
        nuc=NFR_NUC_DIR + "/{celltype}_nucqtl.bed",
        ca=OUTDIR + "/ca_partition_np_peaks/{celltype}_caqtl.bed",
    conda:
        "atac"
    shell:
        """
        mkdir -p {OUTDIR}/upset_figs_f50
        python ../../scripts/bin/upset_qtl.py {input.ca_nfr} {input.ca_nuc} {input.nfr_nuc} {input.nfr_ca} \
                {input.nuc_ca} {input.nuc_nfr} {input.three_beds} \
                {params.celltype} {params.nfr} {params.nuc} {params.ca} {output}
        """
