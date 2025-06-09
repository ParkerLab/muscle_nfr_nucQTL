configfile: "peakcall.yaml"

celltypes=config['celltypes']

rule all:
    input:
        #expand("../../results/results_04262024/macs2_nucleoatac_comparison/nfr/{celltype}_macs2_in_nucleoatac_nfr.bed", celltype=celltypes),
        #expand("../../results/results_04262024/macs2_nucleoatac_comparison/nuc/{celltype}_macs2_in_nucleoatac_nuc.bed", celltype=celltypes),
        #expand("../../results/results_04052024/macs2_nucleoatac_comparison/nuc_narrowpeak/{celltype}_macs2_in_nucleoatac_nuc.bed", celltype=celltypes)
        expand("../../results/results_12032024/summits_distance/{celltype}/{celltype}_nuc_nfr_summits_distance.txt", celltype=celltypes)

# use this rule to compare with nucleoatac inferred nfr regions
rule vs_nucleoatac_nfr:
    input:
        macs2_nfr="../../results/results_04102024/nfrPeaks/nfr_{celltype}.bed",
        nucleoatac_nfr="../../results/results_12192023/nucleoatac/{celltype}/{celltype}.nfr.bed"
    params:
        celltype = "{celltype}",
        #macs2_nfr_bed="../../results/results_03292024/macs2_nfr/{celltype}/atac_macs2_nfr_{celltype}_narrowPeak.bed"
    output:
        dirA="../../results/results_04262024/macs2_nucleoatac_comparison/nfr/{celltype}_macs2_in_nucleoatac_nfr.bed",
        dirB="../../results/results_04262024/macs2_nucleoatac_comparison/nfr/{celltype}_nucleoatac_in_macs2_nfr.bed",
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        # bedtools intersect based on nucleoatac results
        bedtools intersect -a {input.macs2_nfr} -b {input.nucleoatac_nfr} -u > {output.dirA}
        # bedtools intersect based on macs2 results
        bedtools intersect -a {input.nucleoatac_nfr} -b {input.macs2_nfr} -u > {output.dirB}
        """

# use this rule to compare with nucleoatac inferred nuc regions
rule vs_nucleoatac_nuc:
    input:
        #macs2_nuc="../../results/results_04022024/nucPeaks/nuc_{celltype}.bed",
        macs2_nuc="../../results/results_04102024/nucPeaks/nuc_{celltype}.bed",
        nucleoatac_nuc="../../results/results_12192023/nucleoatac/{celltype}/{celltype}.nuc.bed"
    params:
        celltype = "{celltype}",
        #macs2_nuc_bed="../../results/results_03272024/macs2_nuc/atac_macs2_nuc_{celltype}_narrowPeak.bed"
    output:
        dirA="../../results/results_04262024/macs2_nucleoatac_comparison/nuc/{celltype}_macs2_in_nucleoatac_nuc.bed",
        dirB="../../results/results_04262024/macs2_nucleoatac_comparison/nuc/{celltype}_nucleoatac_in_macs2_nuc.bed"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        # bedtools intersect based on nucleoatac results
        bedtools intersect -a {input.macs2_nuc} -b {input.nucleoatac_nuc} -u > {output.dirA}
        # bedtools intersect based on macs2 results
        bedtools intersect -a {input.nucleoatac_nuc} -b {input.macs2_nuc} -u > {output.dirB}
        """

# get the distance between nfr and nuc summits
rule summits_distance:
    input:
        nfr_summits="../../results/results_07112024/summits/{celltype}_nfr_summits.bed",
        nuc_summits="../../results/results_07112024/summits/{celltype}_nuc_summits.bed"
    output:
        "../../results/results_12032024/summits_distance/{celltype}/{celltype}_nuc_nfr_summits_distance.txt"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        bedtools closest -a {input.nuc_summits} -b {input.nfr_summits} -D ref > {output}
        """


    