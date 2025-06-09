configfile: "peakcall.yaml"

celltypes=config['celltypes']
ids=config["IDS"]

import numpy as np
import random

random.seed(42)
np.random.seed(42)


rule all:
    input:
        #expand("../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bed", celltype=celltypes),
        #expand("../../results/results_11212023/macs2/{celltype}/atac_macs2_{celltype}_peaks.broadPeak", celltype=celltypes)
        #expand("../../results/results_07092024/nucleoatac/{celltype}/{celltype}.nucmap_combined.bed.gz", celltype=celltypes),
        #expand("../../results/results_12192023/nucleoatac/{celltype}/{celltype}.nfr.bed", celltype=celltypes)
        #expand("../../data/bams-by-cluster/{celltype}/bychr/atac-{celltype}-fusion_chr{id}.bam", celltype=celltypes, id=ids)
        #expand("../../data/bams-by-cluster/{celltype}/bychr/atac-{celltype}-fusion_chr{id}_sorted.bam", celltype=celltypes, id=ids)
        #expand("../../results/results_07092024/nucleoatac_chr/{celltype}/chr{id}/{celltype}_chr{id}.nucmap_combined.bed.gz", celltype=celltypes, id=ids)
        #expand("../../results/results_02212024/nucleoatac_chr_nfr/{celltype}/chr{id}/{celltype}_chr{id}.nfrpos.bed.gz", celltype=celltypes, id=ids)
        #expand("../../results/results_03292024/macs2_nfr/{celltype}/atac_macs2_nfr_{celltype}_peaks.narrowPeak", celltype=celltypes)
        #expand("../../results/results_03202024/macs2_nfr_100ext/{celltype}/atac_macs2_nfr_{celltype}_peaks.narrowPeak", celltype=celltypes)
        #expand("../../results/results_06282024/macs3_nuc_cutoff0.1_l1_g1/atac_macs2_nuc_{celltype}.narrowPeak", celltype=celltypes)
        #expand("../../results/results_03262024/macs3_cutoff/{celltype}/atac_macs3_nuc_{celltype}_cutoffanalysis", celltype=celltypes)
        #expand("../../results/results_07112024/nucPeaks/nuc_{celltype}.bed", celltype=celltypes),
        #expand("../../results/results_07112024/nfrPeaks/nfr_{celltype}.bed", celltype=celltypes),
        #expand("../../results/results_08042024/macs2/{celltype}/atac_macs2_{celltype}_peaks.narrowPeak", celltype=celltypes)
        expand("../../results/results_10212024/macs2_100ext/{celltype}/atac_macs2_{celltype}_treat_pileup.bw", celltype=celltypes)
        #expand("../../results/results_04082024/macs2_nuc_bybed_summits/{celltype}/atac_macs2_nuc_{celltype}_peaks.narrowPeak", celltype=celltypes)
        #expand("../../results/results_04082024/macs2_nuc_bdg/atac_macs2_nuc_{celltype}_cutoff_analaysis.txt", celltype=celltypes)

rule bamtobed:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bam"
    params:
        celltype = "{celltype}"
    output:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bed"
    envmodules:
        "Bioinformatics",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        bedtools bamtobed -i {input} > {output}
        """

rule macs2:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bed"
    params:
        celltype = "{celltype}"
    output:
        "../../results/results_11212023/macs2/{celltype}/atac_macs2_{celltype}_peaks.broadPeak"
    conda:
        "atac"
    shell:
        """
        macs2 callpeak \
	        -t {input} \
	        --outdir ../../results/results_11212023/macs2/{params.celltype} \
	        --bdg \
	        --SPMR \
	        -f BED \
	        -n atac_macs2_{params.celltype} \
	        -g hs \
	        --shift -100 \
	        --extsize 200 \
	        --seed 762873 \
	        --broad \
	        --keep-dup all \
            --nomodel
        """

rule macs2_narrow:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bed"
    params:
        celltype = "{celltype}"
    output:
        "../../results/results_08042024/macs2/{celltype}/atac_macs2_{celltype}_peaks.narrowPeak"
    conda:
        "atac"
    shell:
        """
        macs2 callpeak \
	        -t {input} \
	        --outdir ../../results/results_08042024/macs2/{params.celltype} \
	        --bdg \
	        --SPMR \
	        -f BED \
	        -n atac_macs2_{params.celltype} \
	        -g hs \
	        --shift -50 \
	        --extsize 100 \
	        --seed 762873 \
	        --keep-dup all \
            --nomodel \
            --call-summits
        """

rule macs2_bw:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bed"
    params:
        celltype = "{celltype}",
        sizes = "/scratch/scjp_root/scjp1/xiaoouw/Muscle_snATAC_nucQTL/data/hg38/hg38.chrom.sizes",
        bdg="../../results/results_10212024/macs2_100ext/{celltype}/atac_macs2_{celltype}_treat_pileup.bdg"
    output:
        "../../results/results_10212024/macs2_100ext/{celltype}/atac_macs2_{celltype}_treat_pileup.bw"
    envmodules:
        "Bioinformatics",
        "UCSC-utilities/4.0.0"
    conda:
        "atac"
    shell:
        """
        macs2 callpeak \
	        -t {input} \
	        --outdir ../../results/results_10212024/macs2_100ext/{params.celltype} \
	        --bdg \
	        --SPMR \
	        -f BED \
	        -n atac_macs2_{params.celltype} \
	        -g hs \
	        --shift -100 \
	        --extsize 200 \
	        --seed 762873 \
	        --keep-dup all \
            --nomodel \
            --call-summits
        bedGraphToBigWig {params.bdg} {params.sizes} {output}
        """



rule macs2_nfr:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr.bam"
    params:
        celltype = "{celltype}",
        bed = "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr.bed"
    output:
        "../../results/results_03292024/macs2_nfr/{celltype}/atac_macs2_nfr_{celltype}_peaks.narrowPeak"
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        # bedtools bamtobed -i {input} > {params.bed}
        macs2 callpeak \
	        -t {params.bed} \
	        --outdir ../../results/results_03292024/macs2_nfr/{params.celltype} \
	        --bdg \
	        --SPMR \
	        -f BED \
	        -n atac_macs2_nfr_{params.celltype} \
	        -g hs \
	        --shift -50 \
	        --extsize 100 \
	        --seed 762873 \
	        --keep-dup all \
            --nomodel
        """

rule macs2_nfr_ext:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr.bam"
    params:
        celltype = "{celltype}",
        bed = "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr.bed"
    output:
        "../../results/results_03202024/macs2_nfr_100ext/{celltype}/atac_macs2_nfr_{celltype}_peaks.narrowPeak"
    conda:
        "atac"
    shell:
        """
        macs2 callpeak \
	        -t {params.bed} \
	        --outdir ../../results/results_03202024/macs2_nfr_100ext/{params.celltype} \
	        --bdg \
	        --SPMR \
	        -f BED \
	        -n atac_macs2_nfr_{params.celltype} \
	        -g hs \
	        --shift -100 \
	        --extsize 200 \
	        --seed 762873 \
	        --keep-dup all \
            --nomodel
        """

rule macs2_nuc:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nuc.bam"
    params:
        celltype = "{celltype}",
        bed = "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nuc.bed",
    output:
        "../../results/results_04082024/macs2_nuc_bybed_summits/{celltype}/atac_macs2_nuc_{celltype}_peaks.narrowPeak"
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        # bedtools bamtobed -i {input} > {params.bed}
        macs2 callpeak \
	        -t {params.bed} \
	        --outdir ../../results/results_04082024/macs2_nuc_bybed_summits/{params.celltype} \
	        --bdg \
	        --SPMR \
	        -f BED \
	        -n atac_macs2_nuc_{params.celltype} \
	        -g hs \
	        --shift -75 \
	        --extsize 150 \
	        --seed 762873 \
	        --keep-dup all \
            --nomodel \
            --call-summits
        """

rule macs2_nuc_bdg:                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_diff.bdg"
    params:
        celltype = "{celltype}",
        prefix= "../../results/results_04082024/macs2_nuc_bdg/atac_macs3_nuc_{celltype}",
        diff_mod = "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_diff_mod.bdg"
    output:
        #"../../results/results_04082024/macs2_nuc/atac_macs2_nuc_{celltype}.narrowPeak"
        #"../../results/results_04082024/macs2_nuc_bdg/atac_macs2_nuc_{celltype}_cutoff_analaysis.txt"
        "../../results/results_06282024/macs3_nuc_cutoff0.1_l1_g1/atac_macs2_nuc_{celltype}.narrowPeak"
    conda:
        "atac"
    shell:
        """
        #macs2 bdgopt -i {input} -m max -p 0 -o {params.diff_mod}
        macs3 bdgpeakcall -i {params.diff_mod} -o {output} -l 1 -g 1 -c 0.1
        """

rule nucleoatac:
    input:
        "../../results/results_11212023/macs2/{celltype}/atac_macs2_{celltype}_peaks.broadPeak"
    params:
        celltype = "{celltype}"
    output:
        "../../results/results_07092024/nucleoatac/{celltype}/{celltype}.nucmap_combined.bed.gz"
    conda:
        "py27"
    shell:
        """
        nucleoatac run --bed ../../results/results_11212023/macs2/{params.celltype}/atac_macs2_{params.celltype}_peaks.broadPeak \
                        --bam  ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion.bam \
                        --fasta ../../data/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta \
                        --out ../../results/results_07092024/nucleoatac/{params.celltype}/{params.celltype} \
                        --cores 10
        """
rule nucleoatac_BED:
    input:
        "../../results/results_12192023/nucleoatac/{celltype}/{celltype}.nfrpos.bed"
    params:
        celltype = "{celltype}"
    output:
        "../../results/results_12192023/nucleoatac/{celltype}/{celltype}.nfr.bed"
    shell:
        """
        awk '{{print $1"\t"$2-100"\t"$3+100}}' \
            ../../results/results_12192023/nucleoatac/{params.celltype}/{params.celltype}.nucmap_combined.bed > \
            ../../results/results_12192023/nucleoatac/{params.celltype}/{params.celltype}.nuc.bed

        awk '{{print $1"\t"$2"\t"$3}}' \
            ../../results/results_12192023/nucleoatac/{params.celltype}/{params.celltype}.nfrpos.bed > \
            ../../results/results_12192023/nucleoatac/{params.celltype}/{params.celltype}.nfr.bed
        """

rule split_bam:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bam"
    params:
        celltype = "{celltype}",
        chr="chr{id}",
    output:
        "../../data/bams-by-cluster/{celltype}/bychr/atac-{celltype}-fusion_chr{id}.bam"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        samtools view -b {input} {params.chr} > {output}
        """

rule sortindex_bam:
    input:
        "../../data/bams-by-cluster/{celltype}/bychr/atac-{celltype}-fusion_chr{id}.bam"
    output:
        "../../data/bams-by-cluster/{celltype}/bychr/atac-{celltype}-fusion_chr{id}_sorted.bam"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        samtools sort {input} > {output}
        samtools index {output}
        """

rule nucleoatac_chr:
    input:
        "../../data/bams-by-cluster/{celltype}/bychr/atac-{celltype}-fusion_chr{id}.bam"
    params:
        celltype = "{celltype}",
        chr="chr{id}",
        #sortbam = "../../data/bams-by-cluster/{celltype}/bychr/atac-{celltype}-fusion_chr{id}_sorted.bam",
        #sortbai = "../../data/bams-by-cluster/{celltype}/bychr/atac-{celltype}-fusion_chr{id}_sorted.bai"
    output:
        "../../results/results_07092024/nucleoatac_chr/{celltype}/chr{id}/{celltype}_chr{id}.nucmap_combined.bed.gz"
    conda:
        "py27"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        samtools index {input}
        nucleoatac run --bed ../../results/results_11212023/macs2/{params.celltype}/atac_macs2_{params.celltype}_peaks.broadPeak \
                        --bam {input} \
                        --fasta ../../data/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta \
                        --out ../../results/results_07092024/nucleoatac_chr/{params.celltype}/{params.chr}/{params.celltype}_{params.chr} \
                        --cores 10
        """ 

rule nucleoatac_nfr_chr:
    input:
        "../../results/results_11212023/macs2/{celltype}/atac_macs2_{celltype}_peaks.broadPeak"
    params:
        celltype = "{celltype}",
        chr="chr{id}",
        run_output = "../../results/results_01232024/nucleoatac_chr/{celltype}/chr{id}/{celltype}_chr{id}"
    output:
        "../../results/results_02212024/nucleoatac_chr_nfr/{celltype}/chr{id}/{celltype}_chr{id}.nfrpos.bed.gz"
    conda:
        "nucleoatac2"
    shell:
        """
        nucleoatac nfr --bed {input} \
                        --bam ../../data/bams-by-cluster/{params.celltype}/bychr/atac-{params.celltype}-fusion_{params.chr}_sorted.bam \
                        --fasta ../../data/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta \
                        --out ../../results/results_02212024/nucleoatac_chr_nfr/{params.celltype}/{params.chr}/{params.celltype}_{params.chr} \
                        --cores 10 \
                        --max_occ 0.2 \
                        --occ_track {params.run_output}.occ.bedgraph.gz \
                        --call {params.run_output}.nucpos.bed.gz
        """  

rule buildpeaks:
    input:
        nfr_summits="../../results/results_03282024/macs2_nfr/{celltype}/atac_macs2_nfr_{celltype}_summits.bed",
        nuc_summits_raw="../../results/results_12192023/nucleoatac/{celltype}/{celltype}.nucmap_combined.bed"
    output:
        nfrPeaks="../../results/results_07112024/nfrPeaks/nfr_{celltype}.bed",
        nucPeaks="../../results/results_07112024/nucPeaks/nuc_{celltype}.bed"
    params:
        nuc_summits_tmp="../../results/results_07112024/summits/{celltype}_nuc_summits_nucleoatac.bed",
        nuc_summits="../../results/results_07112024/summits/{celltype}_nuc_summits.bed",
        nfr_summits="../../results/results_07112024/summits/{celltype}_nfr_summits.bed",
        dir_summits="../../results/results_07112024/summits/"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm"
    shell:
        """
        mkdir -p {params.dir_summits}
        # nfr
        awk '{{print $1"\t"$2-75"\t"$3+75}}' \
            {input.nfr_summits} > {output.nfrPeaks}
        scp {input.nfr_summits} {params.nfr_summits}

        # nuc summits
        awk '{{print $1"\t"$2"\t"$3}}' \
            {input.nuc_summits_raw} > {params.nuc_summits_tmp}

        bedtools subtract -a {params.nuc_summits_tmp} -b {output.nfrPeaks} > {params.nuc_summits}
        
        # nuc
        awk '{{print $1"\t"$2-100"\t"$3+100}}' \
            {params.nuc_summits} > {output.nucPeaks}
        """