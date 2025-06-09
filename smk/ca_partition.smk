configfile: "ca_partition.yaml"

celltypes=config["celltypes"]


rule all:
    input:
        expand("../../results/results_02262025/ca_partition_{celltype}/all_{celltype}.bed", celltype=celltypes),
        expand("../../results/results_02262025/upset_figs/{celltype}_upsetqtl.png", celltype=celltypes),
        #expand("../../results/results_02262025/ca_partition_figs/ca_partition_{celltype}_venn.png", celltype=celltypes)
        expand("../../results/results_04012025/summits_coverage/summits_in_nuc_{celltype}.bed", celltype=celltypes),
        expand("../../results/results_04012025/summits_coverage/summits_in_nfr_{celltype}.bed", celltype=celltypes),
        expand("../../results/results_04132025/downsample200k_{celltype}/all_{celltype}.bed", celltype=celltypes),
        expand("../../results/results_04132025/upset_figs/{celltype}_upsetqtl.png", celltype=celltypes),
rule intersect:
    input:
        nfr="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_nfrqtl.bed",
        nuc="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_nucqtl.bed",
        ca="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_caqtl.bed",
    output:
        ca_nfr="../../results/results_02262025/ca_partition_{celltype}/ca_nfr_{celltype}.bed",
        ca_nuc="../../results/results_02262025/ca_partition_{celltype}/ca_nuc_{celltype}.bed",
        nfr_nuc="../../results/results_02262025/ca_partition_{celltype}/nfr_nuc_{celltype}.bed",
        nfr_ca="../../results/results_02262025/ca_partition_{celltype}/nfr_ca_{celltype}.bed",
        nuc_ca="../../results/results_02262025/ca_partition_{celltype}/nuc_ca_{celltype}.bed",
        nuc_nfr="../../results/results_02262025/ca_partition_{celltype}/nuc_nfr_{celltype}.bed",
        three_beds="../../results/results_02262025/ca_partition_{celltype}/all_{celltype}.bed"
    params:
        temp="../../results/results_02262025/ca_partition_{celltype}/nfr_nuc_precise_{celltype}.bed"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm",
    shell:
        """
        #bedtools intersect -wa -a {input.ca} -b {input.nfr} -f 0.50 -r > {output.ca_nfr}
        #bedtools intersect -wa -a {input.ca} -b {input.nuc} -f 0.50 -r > {output.ca_nuc}
        #bedtools intersect -wa -a {input.nfr} -b {input.nuc} -f 0.50 -r > {output.nfr_nuc}
        #bedtools intersect -wa -a {input.nfr} -b {input.ca} -f 0.50 -r > {output.nfr_ca}
        #bedtools intersect -wa -a {input.nuc} -b {input.ca} -f 0.50 -r > {output.nuc_ca}
        #bedtools intersect -wa -a {input.nuc} -b {input.nfr} -f 0.50 -r > {output.nuc_nfr}
        #bedtools intersect -a {input.nfr} -b {input.nuc} -f 0.50 -r > {params.temp}
        #bedtools intersect -wa -a {input.ca} -b {params.temp} -f 0.50 -r > {output.three_beds}
        bedtools intersect -wa -a {input.ca} -b {input.nfr} > {output.ca_nfr}
        bedtools intersect -wa -a {input.ca} -b {input.nuc} > {output.ca_nuc}
        bedtools intersect -wa -a {input.nfr} -b {input.nuc} > {output.nfr_nuc}
        bedtools intersect -wa -a {input.nfr} -b {input.ca} > {output.nfr_ca}
        bedtools intersect -wa -a {input.nuc} -b {input.ca} > {output.nuc_ca}
        bedtools intersect -wa -a {input.nuc} -b {input.nfr} > {output.nuc_nfr}
        bedtools intersect -a {input.nfr} -b {input.nuc} > {params.temp}
        bedtools intersect -wa -a {input.ca} -b {params.temp} > {output.three_beds}
        """

rule intersect_summits:
    input:
        macs2_summits="../../results/results_08042024/macs2/{celltype}/atac_macs2_{celltype}_summits.bed",
        nfrpeaks="../../results/results_07112024/nfrPeaks/nfr_{celltype}.bed",
        nucpeaks="../../results/results_07112024/nucPeaks/nuc_{celltype}.bed"
    output:
        summits_in_nfr="../../results/results_04012025/summits_coverage/summits_in_nfr_{celltype}.bed",
        summits_in_nuc="../../results/results_04012025/summits_coverage/summits_in_nuc_{celltype}.bed",
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm",
    shell:
        """
        bedtools intersect -wa -a {input.macs2_summits} -b {input.nfrpeaks} > {output.summits_in_nfr}
        bedtools intersect -wa -a {input.macs2_summits} -b {input.nucpeaks} > {output.summits_in_nuc}
        """

rule venn:
    input:
        ca_nfr="../../results/results_02262025/ca_partition_{celltype}/ca_nfr_{celltype}.bed",
        ca_nuc="../../results/results_02262025/ca_partition_{celltype}/ca_nuc_{celltype}.bed"
    output:
        "../../results/results_02262025/ca_partition_figs/ca_partition_{celltype}_venn.png"
    params:
        celltype="{celltype}",
        caqtl="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_caqtl.bed"
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/ca_venn.py {input.ca_nfr} {input.ca_nuc} {output} {params.celltype} {params.caqtl}
        """

rule upset:
    input:
        ca_nfr="../../results/results_02262025/ca_partition_{celltype}/ca_nfr_{celltype}.bed",
        ca_nuc="../../results/results_02262025/ca_partition_{celltype}/ca_nuc_{celltype}.bed",
        nfr_nuc="../../results/results_02262025/ca_partition_{celltype}/nfr_nuc_{celltype}.bed",
        nfr_ca="../../results/results_02262025/ca_partition_{celltype}/nfr_ca_{celltype}.bed",
        nuc_ca="../../results/results_02262025/ca_partition_{celltype}/nuc_ca_{celltype}.bed",
        nuc_nfr="../../results/results_02262025/ca_partition_{celltype}/nuc_nfr_{celltype}.bed",
        three_beds="../../results/results_02262025/ca_partition_{celltype}/all_{celltype}.bed"
    output:
        "../../results/results_02262025/upset_figs/{celltype}_upsetqtl.png"
    params:
        celltype="{celltype}",
        nfr="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_nfrqtl.bed",
        nuc="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_nucqtl.bed",
        ca="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_caqtl.bed",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/upset_qtl.py {input.ca_nfr} {input.ca_nuc} {input.nfr_nuc} {input.nfr_ca} \
                {input.nuc_ca} {input.nuc_nfr} {input.three_beds} \
                {params.celltype} {params.nfr} {params.nuc} {params.ca} {output}
        """

rule intersect_downsampling:
    input:
        nfr="../../results/results_09132024/{celltype}_nfrqtl_200k.bed",
        nuc="../../results/results_09132024/{celltype}_nucqtl_200k.bed",
        ca="../../results/results_09132024/{celltype}_caqtl_200k.bed",
    output:
        ca_nfr="../../results/results_04132025/downsample200k_{celltype}/ca_nfr_{celltype}.bed",
        ca_nuc="../../results/results_04132025/downsample200k_{celltype}/ca_nuc_{celltype}.bed",
        nfr_nuc="../../results/results_04132025/downsample200k_{celltype}/nfr_nuc_{celltype}.bed",
        nfr_ca="../../results/results_04132025/downsample200k_{celltype}/nfr_ca_{celltype}.bed",
        nuc_ca="../../results/results_04132025/downsample200k_{celltype}/nuc_ca_{celltype}.bed",
        nuc_nfr="../../results/results_04132025/downsample200k_{celltype}/nuc_nfr_{celltype}.bed",
        three_beds="../../results/results_04132025/downsample200k_{celltype}/all_{celltype}.bed"
    params:
        temp="../../results/results_04132025/downsample200k_{celltype}/nfr_nuc_precise_{celltype}.bed"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm",
    shell:
        """
        bedtools intersect -wa -a {input.ca} -b {input.nfr} -f 0.50 -r > {output.ca_nfr}
        bedtools intersect -wa -a {input.ca} -b {input.nuc} -f 0.50 -r > {output.ca_nuc}
        bedtools intersect -wa -a {input.nfr} -b {input.nuc} -f 0.50 -r > {output.nfr_nuc}
        bedtools intersect -wa -a {input.nfr} -b {input.ca} -f 0.50 -r > {output.nfr_ca}
        bedtools intersect -wa -a {input.nuc} -b {input.ca} -f 0.50 -r > {output.nuc_ca}
        bedtools intersect -wa -a {input.nuc} -b {input.nfr} -f 0.50 -r > {output.nuc_nfr}
        bedtools intersect -a {input.nfr} -b {input.nuc} -f 0.50 -r > {params.temp}
        bedtools intersect -wa -a {input.ca} -b {params.temp} -f 0.50 -r > {output.three_beds}
        """

rule upset_downsampling:
    input:
        ca_nfr="../../results/results_04132025/downsample200k_{celltype}/ca_nfr_{celltype}.bed",
        ca_nuc="../../results/results_04132025/downsample200k_{celltype}/ca_nuc_{celltype}.bed",
        nfr_nuc="../../results/results_04132025/downsample200k_{celltype}/nfr_nuc_{celltype}.bed",
        nfr_ca="../../results/results_04132025/downsample200k_{celltype}/nfr_ca_{celltype}.bed",
        nuc_ca="../../results/results_04132025/downsample200k_{celltype}/nuc_ca_{celltype}.bed",
        nuc_nfr="../../results/results_04132025/downsample200k_{celltype}/nuc_nfr_{celltype}.bed",
        three_beds="../../results/results_04132025/downsample200k_{celltype}/all_{celltype}.bed"
    output:
        "../../results/results_04132025/upset_figs/{celltype}_upsetqtl.png"
    params:
        celltype="{celltype}",
        nfr="../../results/results_09132024/{celltype}_nfrqtl_200k.bed",
        nuc="../../results/results_09132024/{celltype}_nucqtl_200k.bed",
        ca="../../results/results_09132024/{celltype}_caqtl_200k.bed",
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/upset_qtl.py {input.ca_nfr} {input.ca_nuc} {input.nfr_nuc} {input.nfr_ca} \
                {input.nuc_ca} {input.nuc_nfr} {input.three_beds} \
                {params.celltype} {params.nfr} {params.nuc} {params.ca} {output}
        """