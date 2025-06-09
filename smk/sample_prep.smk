configfile: "sample_prep.yaml"

celltypes=config['celltypes']
regions=config['regions']

# get the samples
import pandas as pd

def get_samples(sample_file):
    sample_df = pd.read_csv(sample_file, sep = '\t')
    # get ids
    samples = sample_df['SAMPLE']+'--'+sample_df['batch']
    # get unique ids
    samples = list(samples.drop_duplicates())
    # remove '-NM'
    return([s.strip('-NM') for s in samples])

samples = get_samples("../..data/sample_info/sample_level_covariates_atac.tsv")

rule all:
    input:
        #expand("../../data/bams-by-cluster-sample/{celltype}/chr22/sizes/atac-{celltype}@fusion@{sample}_chr22_sizes.txt", celltype=celltypes,sample=samples)
        #expand("../../data/bams-by-cluster-sample/{celltype}/chr22/thresholds/atac-{celltype}@fusion@{sample}_chr22_thresholds.txt", celltype=celltypes,sample=samples)
        #expand("../../data/bams-by-cluster-sample/{celltype}/nuc/atac-{celltype}@fusion@{sample}_nuc.bam", celltype=celltypes,sample=samples),
        #expand("../../data/bams-by-cluster-sample/{celltype}/nfr/atac-{celltype}@fusion@{sample}_nfr.bam", celltype=celltypes,sample=samples)
        #expand("../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr.bam",celltype=celltypes),
        #expand("../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nuc.bam",celltype=celltypes),
        #expand("../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nuc_sorted.bam",celltype=celltypes),
        #expand("../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr_sorted.bam",celltype=celltypes)
        #expand("../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_{region}_sorted.bed", celltype=celltypes, region=regions)
        #expand("../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_{region}.bw", celltype=celltypes, region=regions)
        #expand("../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bw", celltype=celltypes)
        #expand("../../data/sample_info/{celltype}_sample_level_covariates_atac.tsv", celltype=celltypes)
        #expand("../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_diff.bw", celltype=celltypes)
        #expand("../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_diff.bdg", celltype=celltypes)
        #expand("../../data/bams-by-cluster-sample/{celltype}/{region}_bw/atac-{celltype}@fusion@{sample}_{region}.bw", sample=samples, celltype=celltypes, region=regions),
        #expand("../../data/bams-by-cluster-sample/{celltype}/diff_bw_mod/atac-{celltype}@fusion@{sample}_diff.bw", sample=samples, celltype=celltypes)
        #expand("../../data/bams-by-cluster-sample/{celltype}/full_bw/atac-{celltype}@fusion@{sample}.bw", sample=samples, celltype=celltypes)
        expand("../../data/bams-by-cluster-sample/{celltype}/diff_bw_mod_smooth_25bp/atac-{celltype}@fusion@{sample}_diff.bw", sample=samples, celltype=celltypes)


rule split_chr22:
    input:
        "../../data/bams-by-cluster-sample/{celltype}/atac-{celltype}@fusion@{sample}.bam"
    params:
        sample = "{sample}",
        celltype = "{celltype}"
    output:
        "../../data/bams-by-cluster-sample/{celltype}/chr22/sizes/atac-{celltype}@fusion@{sample}_chr22_sizes.txt"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        samtools view -bo ../../data/bams-by-cluster-sample/{params.celltype}/chr22/atac-{params.celltype}@fusion@{params.sample}_chr22.bam {input} chr22
        samtools index ../../data/bams-by-cluster-sample/{params.celltype}/chr22/atac-{params.celltype}@fusion@{params.sample}_chr22.bam
        samtools view -f66 ../../data/bams-by-cluster-sample/{params.celltype}/chr22/atac-{params.celltype}@fusion@{params.sample}_chr22.bam \
          | cut -f9 | awk '{{print sqrt($0^2)}}' > \
          ../../data/bams-by-cluster-sample/{params.celltype}/chr22/sizes/atac-{params.celltype}@fusion@{params.sample}_chr22_sizes.txt
        """

rule thresholding:
    input:
        "../../data/bams-by-cluster-sample/{celltype}/chr22/sizes/atac-{celltype}@fusion@{sample}_chr22_sizes.txt"
    params:
        sample = "{sample}",
        celltype = "{celltype}"
    output:
        "../../data/bams-by-cluster-sample/{celltype}/chr22/thresholds/atac-{celltype}@fusion@{sample}_chr22_thresholds.txt"
    conda:
        "atac"
    shell:
        """
        python ../bin/multiotsu_threshold.py \
        ../../data/bams-by-cluster-sample/{params.celltype}/chr22/sizes/atac-{params.celltype}@fusion@{params.sample}_chr22_sizes.txt \
        ../../data/bams-by-cluster-sample/{params.celltype}/chr22/thresholds/atac-{params.celltype}@fusion@{params.sample}_chr22_thresholds.txt
        """

rule split_bam:
    input:
        "../../data/bams-by-cluster-sample/{celltype}/chr22/thresholds/atac-{celltype}@fusion@{sample}_chr22_thresholds.txt"
    params:
        sample = "{sample}",
        celltype = "{celltype}"
    output:
        "../../data/bams-by-cluster-sample/{celltype}/nuc/atac-{celltype}@fusion@{sample}_nuc.bam",
        "../../data/bams-by-cluster-sample/{celltype}/nfr/atac-{celltype}@fusion@{sample}_nfr.bam"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        t=$(head -n 1 {input})
        ## filter - NFR
        samtools view -h ../../data/bams-by-cluster-sample/{params.celltype}/atac-{params.celltype}@fusion@{params.sample}.bam \
        | awk -v r=$t '(substr($0,1,1)=="@" || ($9>= -1*r && $9<=r))' | \
        samtools view -b > ../../data/bams-by-cluster-sample/{params.celltype}/nfr/atac-{params.celltype}@fusion@{params.sample}_nfr.bam
        ## index
        samtools index ../../data/bams-by-cluster-sample/{params.celltype}/nfr/atac-{params.celltype}@fusion@{params.sample}_nfr.bam

        ## filter - Nucleosomal
        samtools view -h ../../data/bams-by-cluster-sample/{params.celltype}/atac-{params.celltype}@fusion@{params.sample}.bam \
        | awk -v r=$t '(substr($0,1,1)=="@" || ($9<= -1*r || $9>=r))' | \
        samtools view -b > ../../data/bams-by-cluster-sample/{params.celltype}/nuc/atac-{params.celltype}@fusion@{params.sample}_nuc.bam
        ## index
        samtools index ../../data/bams-by-cluster-sample/{params.celltype}/nuc/atac-{params.celltype}@fusion@{params.sample}_nuc.bam
        """

rule merge_bams:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bam"
    params:
        celltype = "{celltype}"
    output:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr.bam",
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nuc.bam"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        # nfr
        samtools merge ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion_nfr.bam \
        ../../data/bams-by-cluster-sample/{params.celltype}/nfr/*.bam

        # nuc
        samtools merge ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion_nuc.bam \
        ../../data/bams-by-cluster-sample/{params.celltype}/nuc/*.bam

        # index
        samtools index ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion_nfr.bam
        samtools index ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion_nuc.bam
        """

rule sort_bam_nfr:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr.bam"
    params:
        celltype = "{celltype}"
    output:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr_sorted.bam"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        # nfr
        samtools sort -l 1 -n -m 8G -@ 18 ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion_nfr.bam \
        -T ../../data/bams-by-cluster/{params.celltype}/{params.celltype}_sort_tmp_nfr \
        -o ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion_nfr_sorted.bam
        """

rule sort_bam_nuc:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nuc.bam"
    params:
        celltype = "{celltype}"
    output:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nuc_sorted.bam"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        # nuc
        samtools sort -l 1 -n -m 8G -@ 18 ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion_nuc.bam \
        -T ../../data/bams-by-cluster/{params.celltype}/{params.celltype}_sort_tmp_nuc \
        -o ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion_nuc_sorted.bam
        """
rule bam_to_bedpe:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_{region}_sorted.bam"
    params:
        celltype = "{celltype}",
        region = "{region}"
    output:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_{region}_sorted.bed"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n",
        "bedtools2/2.30.0-svcfwbm",
        "UCSC-utilities/4.0.0"
    shell:
        """
        bedtools bamtobed -i ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion_{params.region}_sorted.bam \
        -bedpe > ../../data/bams-by-cluster/{params.celltype}/atac-{params.celltype}-fusion_{params.region}_sorted.bed
        """

rule bamCoverage:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_{region}.bam"
    params:
        celltype = "{celltype}",
        region = "{region}",
        blacklist = "../..data/hg38/ENCFF356LFX.bed.gz"
    output:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_{region}.bw"
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        bamCoverage -b {input} \
            -o {output} \
            --normalizeUsing BPM \
            --numberOfProcessors 8 \
            --blackListFileName {params.blacklist} \
            --extendReads \
            --binSize 10 \
            --smoothLength 30
        """

rule bamCoverage_full:
    input:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bam"
    params:
        celltype = "{celltype}",
        blacklist = "../..data/hg38/ENCFF356LFX.bed.gz"
    output:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion.bw"
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        bamCoverage -b {input} \
            -o {output} \
            --normalizeUsing BPM \
            --numberOfProcessors 8 \
            --blackListFileName {params.blacklist} \
            --extendReads \
            --binSize 10 \
            --smoothLength 30
        """

rule sample_bamCoverage:
    input:
        "../../data/bams-by-cluster-sample/{celltype}/{region}/atac-{celltype}@fusion@{sample}_{region}.bam",
    params:
        celltype = "{celltype}",
        region = "{region}",
        blacklist = "../..data/hg38/ENCFF356LFX.bed.gz",
        sample = "{sample}",
    output:
        "../../data/bams-by-cluster-sample/{celltype}/{region}_bw/atac-{celltype}@fusion@{sample}_{region}.bw",
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        bamCoverage -b {input} \
            -o {output} \
            --normalizeUsing BPM \
            --numberOfProcessors 1 \
            --blackListFileName {params.blacklist} \
            --extendReads \
            --binSize 50 \
            --smoothLength 150
        """

rule sample_bamCoverage_full:
    input:
        "../../data/bams-by-cluster-sample/{celltype}/atac-{celltype}@fusion@{sample}.bam",
    params:
        celltype = "{celltype}",
        blacklist = "../..data/hg38/ENCFF356LFX.bed.gz",
        sample = "{sample}",
    output:
        "../../data/bams-by-cluster-sample/{celltype}/full_bw/atac-{celltype}@fusion@{sample}.bw",
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        bamCoverage -b {input} \
            -o {output} \
            --normalizeUsing BPM \
            --numberOfProcessors 1 \
            --blackListFileName {params.blacklist} \
            --extendReads \
            --binSize 10 \
            --smoothLength 30
        """

rule sample_diff_bw:
    input:
        nfr_bw="../../data/bams-by-cluster-sample/{celltype}/nfr_bw/atac-{celltype}@fusion@{sample}_nfr.bw",
        nuc_bw="../../data/bams-by-cluster-sample/{celltype}/nuc_bw/atac-{celltype}@fusion@{sample}_nuc.bw"
    params:
        celltype= "{celltype}",
        sample = "{sample}",
        blacklist = "../..data/hg38/ENCFF356LFX.bed.gz",
        raw_bdg = "../../data/bams-by-cluster-sample/{celltype}/diff_bdg_raw/atac-{celltype}@fusion@{sample}_diff.bdg",
        modifid_bdg = "../../data/bams-by-cluster-sample/{celltype}/diff_bdg_mod/atac-{celltype}@fusion@{sample}_diff.bdg",
        sorted_bdg = "../../data/bams-by-cluster-sample/{celltype}/diff_bdg_mod/atac-{celltype}@fusion@{sample}_diff_sorted.bdg",
        sizes = "../../data/hg38/hg38.chrom.sizes"
    output:
        "../../data/bams-by-cluster-sample/{celltype}/diff_bw_mod/atac-{celltype}@fusion@{sample}_diff.bw"
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n",
        "UCSC-utilities/4.0.0"
    shell:
        """
        mkdir -p ../../data/bams-by-cluster-sample/{params.celltype}/diff_bdg_raw/
        mkdir -p ../../data/bams-by-cluster-sample/{params.celltype}/diff_bdg_mod/

        bigwigCompare -b1 {input.nuc_bw} -b2 {input.nfr_bw} \
            --operation subtract \
            --skipNAs \
            --blackListFileName {params.blacklist} \
            -p 1 \
            --verbose \
            --outFileName {params.raw_bdg} \
            --outFileFormat bedgraph
        
        macs2 bdgopt -i {params.raw_bdg} -m max -p 0 -o {params.modifid_bdg}
        # sort bdg
        #sort -k1,1 -k2,2n {params.modifid_bdg} > {params.sorted_bdg}
        {{ head -n 1 {params.modifid_bdg}; tail -n +2 {params.modifid_bdg} | sort -k1,1 -k2,2n; }} > {params.sorted_bdg}



        bedGraphToBigWig {params.sorted_bdg} {params.sizes} {output}
        """
    
rule get_chr22_median:
    input:
        "../../data/sample_info/sample_level_covariates_atac.tsv"
    params:
        celltype = "{celltype}"
    output:
        "../../data/sample_info/{celltype}_sample_level_covariates_atac.tsv"
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/get_chr22_median.py {params.celltype} {input} {output}
        """

rule bw_subtraction:
    input:
        nfr="../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr.bw",
        nuc="../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nuc.bw"
    params:
        celltype = "{celltype}",
        blacklist = "../..data/hg38/ENCFF356LFX.bed.gz"
    output:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_diff.bw"
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        bigwigCompare -b1 {input.nuc} -b2 {input.nfr} \
            --operation subtract \
            --skipNAs \
            --blackListFileName {params.blacklist} \
            -p 8 \
            --verbose \
            --outFileName {output} \
            --outFileFormat bigwig
        """

rule bw_subtraction_bdg:
    input:
        nfr="../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nfr.bw",
        nuc="../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_nuc.bw"
    params:
        celltype = "{celltype}",
        blacklist = "../..data/hg38/ENCFF356LFX.bed.gz"
    output:
        "../../data/bams-by-cluster/{celltype}/atac-{celltype}-fusion_diff.bdg"
    conda:
        "atac"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        bigwigCompare -b1 {input.nuc} -b2 {input.nfr} \
            --operation subtract \
            --skipNAs \
            --blackListFileName {params.blacklist} \
            -p 8 \
            --verbose \
            --outFileName {output} \
            --outFileFormat bedgraph
        """

rule sample_nuc_bw_smooth:
    input:
        "../../data/bams-by-cluster-sample/{celltype}/diff_bw_mod/atac-{celltype}@fusion@{sample}_diff_sorted.bdg"
    output:
        "../../data/bams-by-cluster-sample/{celltype}/diff_bw_mod_smooth_25bp/atac-{celltype}@fusion@{sample}_diff.bdg"
    conda:
        "atac"
    shell:
        """
        python ../../scripts/bin/smooth_bdg.py {input} {output}
        """
rule smooth_bdg_sort_to_bw:
    input:
        "../../data/bams-by-cluster-sample/{celltype}/diff_bw_mod_smooth_25bp/atac-{celltype}@fusion@{sample}_diff.bdg"
    output:
        "../../data/bams-by-cluster-sample/{celltype}/diff_bw_mod_smooth_25bp/atac-{celltype}@fusion@{sample}_diff.bw"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n",
        "bedtools2/2.30.0-svcfwbm",
        "UCSC-utilities/4.0.0"
    params:
        sort_bdg="../../data/bams-by-cluster-sample/{celltype}/diff_bw_mod_smooth_25bp/atac-{celltype}@fusion@{sample}_diff_sort.bdg",
        sizes="../../data/hg38/hg38.chrom.sizes"
    shell:
        """
        # sort bdg
        {{ head -n 1 {input}; tail -n +2 {input} | sort -k1,1 -k2,2n; }} > {params.sort_bdg}
        # convert to bw
        bedGraphToBigWig {params.sort_bdg} {params.sizes} {output}
        """