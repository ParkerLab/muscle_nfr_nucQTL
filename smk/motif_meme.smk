configfile: "motif_meme.yaml"

celltypes=config["celltypes"]
regions=config["regions"]
ids=config["IDS"]

rule all:
    input:
        #expand("../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_nfrqtl.bed", celltype=celltypes),
        #expand("../../results/results_10102024/motif_meme_qtl/seqs/{celltype}_{region}qtl.fasta", celltype=celltypes, region=regions),
        #expand("../../results/results_12032024/motif_meme_qtl/sea_{celltype}_{region}/sea.tsv", celltype=celltypes, region=regions),
        #expand("../../results/results_02182025/motif_meme_coloc/seqs/{celltype}_nfrqtl_no_coloc.fasta", celltype=celltypes),
        #expand("../../results/results_02182025/motif_meme_coloc/sea_{celltype}/sea.tsv", celltype=celltypes),
        #expand("../../results/results_11252024/motif_meme_peak/{celltype}_caPeak.fasta", celltype=celltypes),
        #expand("../../results/results_12032024/motif_meme_peak/sea_{celltype}_{region}Peak/sea.tsv", celltype=celltypes, region=regions),
        #expand("../../results/results_01022025/motif_meme_snpaware/{celltype}/chr{id}_close.fasta", celltype=celltypes, id=ids),
        #expand("../../results/results_01032025/motif_meme_snpaware/seqs/{celltype}_nfrqtl_close.fasta", celltype=celltypes),
        expand("../../results/results_03062025/motif_meme_snpaware_ctr/seqs/{celltype}_nfrqtl_close.fasta", celltype=celltypes),
        expand("../../results/results_03062025/motif_meme_snpaware_ctr/seqs/{celltype}_nfrqtl_open.fasta", celltype=celltypes),
        #expand('../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/chr22_open.fasta', celltype=celltypes),
        #expand('../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/chr22_close.fasta', celltype=celltypes),
        #expand("../../results/results_01032025/motif_meme_snpaware/sea_{celltype}/sea.tsv", celltype=celltypes)
        expand("../../results/results_03062025/motif_meme_snpaware_ctr/sea_{celltype}/sea.tsv", celltype=celltypes),
        #expand("../../results/results_04232025/motif_meme_cit/sea_{celltype}/sea.tsv", celltype=celltypes),
        expand('../../results/results_04232025/motif_meme_cit/sea_merged/sea.tsv'),


rule formatnfr:
    input:
        nfrcoloc="../../results/results_02182025/motif_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_nuc.homerpeak",
        nfrno="../../results/results_02182025/motif_coloc/{celltype}_nfrqtl_coloc/{celltype}_nfr_coloc_nuc_bkg.homerpeak"
    output:
        nfrcoloc="../../results/results_02182025/motif_meme_coloc/peaks/{celltype}_nfrqtl_coloc.bed",
        nfrno="../../results/results_02182025/motif_meme_coloc/peaks/{celltype}_nfrqtl_no_coloc.bed"
    shell:
        """
        awk 'OFS="\t" {{print $2, $3, $4}}' {input.nfrcoloc} > {output.nfrcoloc}
        awk 'OFS="\t" {{print $2, $3, $4}}' {input.nfrno} > {output.nfrno}
        """

rule formatqtl:
    input:
        nfr="../../results/results_05102024/QTL_opt_results/{celltype}_nfr/merged_QTLresults_fdrcorr.csv",
        nuc="../../results/results_07162024/QTL_opt_results/{celltype}_nuc/merged_QTLresults_fdrcorr.csv",
        ca='../../results/results_08312024/QTL_opt_results/{celltype}_ca_merged_QTLresults_fdrcorr.csv',
    output:
        nfr="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_nfrqtl.bed",
        nuc="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_nucqtl.bed",
        ca="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_caqtl.bed",
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/qtl_csv_to_bed.R --csv {input.nfr} --out {output.nfr}
        Rscript ../../scripts/bin/qtl_csv_to_bed.R --csv {input.nuc} --out {output.nuc}
        Rscript ../../scripts/bin/qtl_csv_to_bed.R --csv {input.ca} --out {output.ca}
        """

rule getfasta_qtl:
    input:
        hg38="../../data/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta",
        qtlbed="../../results/results_10102024/motif_meme_qtl/peaks/{celltype}_{region}qtl.bed",
    output:
        "../../results/results_10102024/motif_meme_qtl/seqs/{celltype}_{region}qtl.fasta",
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm",
        "samtools/1.13-fwwss5n"
    shell:
        """
        bedtools getfasta -fi {input.hg38} -bed {input.qtlbed} -fo {output}
        samtools faidx {output}
        """

rule getfasta_peak:
    input:
        nfr="../../results/results_07112024/nfrPeaks/nfr_{celltype}.bed",
        nuc="../../results/results_07112024/nucPeaks/nuc_{celltype}.bed",
        ca="../../results/results_08292024/caPeaks/ca_{celltype}_75.bed"
    params:
        hg38="../../data/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
    output:
        nfr="../../results/results_11252024/motif_meme_peak/{celltype}_nfrPeak.fasta",
        nuc="../../results/results_11252024/motif_meme_peak/{celltype}_nucPeak.fasta",
        ca="../../results/results_11252024/motif_meme_peak/{celltype}_caPeak.fasta"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm",
        "samtools/1.13-fwwss5n"
    shell:
        """
        bedtools getfasta -fi {params.hg38} -bed {input.nfr} -fo {output.nfr}
        samtools faidx {output.nfr}
        bedtools getfasta -fi {params.hg38} -bed {input.nuc} -fo {output.nuc}
        samtools faidx {output.nuc}
        bedtools getfasta -fi {params.hg38} -bed {input.ca} -fo {output.ca}
        samtools faidx {output.ca}
        """


rule getfasta_coloc:
    input:
        hg38="../../data/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta",
        nfrcoloc="../../results/results_02182025/motif_meme_coloc/peaks/{celltype}_nfrqtl_coloc.bed",
        nfrno="../../results/results_02182025/motif_meme_coloc/peaks/{celltype}_nfrqtl_no_coloc.bed",
    output:
        nfrcoloc="../../results/results_02182025/motif_meme_coloc/seqs/{celltype}_nfrqtl_coloc.fasta",
        nfrno="../../results/results_02182025/motif_meme_coloc/seqs/{celltype}_nfrqtl_no_coloc.fasta",
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm",
        "samtools/1.13-fwwss5n"
    shell:
        """
        bedtools getfasta -fi {input.hg38} -bed {input.nfrcoloc} -fo {output.nfrcoloc}
        samtools faidx {output.nfrcoloc}
        bedtools getfasta -fi {input.hg38} -bed {input.nfrno} -fo {output.nfrno}
        samtools faidx {output.nfrno}
        """

rule getfasta_cit:
    input:
        hg38="../../data/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta",
        nfrcit="../../results/results_07212024/{celltype}_cit_nfr_to_nuc.bed",
    output:
        nfrcit="../../results/results_04232025/motif_meme_cit/seqs/{celltype}_cit.fasta"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm",
        "samtools/1.13-fwwss5n"
    shell:
        """
        bedtools getfasta -fi {input.hg38} -bed {input.nfrcit} -fo {output.nfrcit}
        samtools faidx {output.nfrcit}
        """

rule getfasta_snpaware:
    input:
        "../../results/results_11182024/{celltype}_snps_inpeak_info/{celltype}_chr{id}_nfr_cs_pip.tsv"
    params:
        var_info="../../data/sample_info/variant_info.tsv",
        nfrqtl_fasta="../../results/results_10102024/motif_meme_coloc/seqs/{celltype}_nfrqtl_coloc.fasta",
    output:
        open="../../results/results_01022025/motif_meme_snpaware/{celltype}/chr{id}_open.fasta",
        close="../../results/results_01022025/motif_meme_snpaware/{celltype}/chr{id}_close.fasta"
    conda:
        "R423"
    shell:
        """
        Rscript ../../scripts/bin/get_snpaware_fasta.R \
            --inpeak_dir {input} \
            --var_info {params.var_info} \
            --nfrqtl_fasta {params.nfrqtl_fasta} \
            --output_open {output.open} \
            --output_close {output.close}
        """

rule getfasta_snpaware_ctr:
    input:
        "../../results/results_11182024/{celltype}_snps_inpeak_info/{celltype}_chr{id}_nfr_cs_pip.tsv"
    params:
        var_info="../../data/sample_info/variant_info.tsv",
        hg38_fasta="../../data/hg38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta",
        part1_bed="../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/chr{id}_tmp_part1.bed",
        part2_bed="../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/chr{id}_tmp_part2.bed",
        part1_fasta="../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/chr{id}_tmp_part1.fasta",
        part2_fasta="../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/chr{id}_tmp_part2.fasta",
    output:
        open="../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/chr{id}_open.fasta",
        close="../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/chr{id}_close.fasta"
    conda:
        "R423"
    envmodules:
        "Bioinformatics",
        "gcc/10.3.0-k2osx5y",
        "bedtools2/2.30.0-svcfwbm",
        "samtools/1.13-fwwss5n"
    shell:
        """
        Rscript ../../scripts/bin/get_snpaware_fasta_ctr.R \
            --inpeak_dir {input} \
            --var_info {params.var_info} \
            --hg38_fasta {params.hg38_fasta} \
            --part1_bed {params.part1_bed} \
            --part2_bed {params.part2_bed} \
            --part1_fasta {params.part1_fasta} \
            --part2_fasta {params.part2_fasta} \
            --output_open {output.open} \
            --output_close {output.close}
        rm -r {params.part1_bed}
        rm -r {params.part2_bed}
        rm -r {params.part1_fasta}
        rm -r {params.part2_fasta}
        """

rule merge_snpaware_fasta:
    input:
        open="../../results/results_01022025/motif_meme_snpaware/{celltype}/chr22_open.fasta",
        close="../../results/results_01022025/motif_meme_snpaware/{celltype}/chr22_close.fasta"
    params:
        open_dir="../../results/results_01022025/motif_meme_snpaware/{celltype}/",
        close_dir="../../results/results_01022025/motif_meme_snpaware/{celltype}/"
    output:
        open_merge="../../results/results_01032025/motif_meme_snpaware/seqs/{celltype}_nfrqtl_open.fasta",
        close_merge="../../results/results_01032025/motif_meme_snpaware/seqs/{celltype}_nfrqtl_close.fasta"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        cat {params.open_dir}chr*_open.fasta > {output.open_merge}
        cat {params.close_dir}chr*_close.fasta > {output.close_merge}
        samtools faidx {output.open_merge}
        samtools faidx {output.close_merge}
        """

rule merge_snpaware_fasta_ctr:
    input:
        open="../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/chr22_open.fasta",
        close="../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/chr22_close.fasta"
    params:
        open_dir="../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/",
        close_dir="../../results/results_03062025/motif_meme_snpaware_ctr/{celltype}/"
    output:
        open_merge="../../results/results_03062025/motif_meme_snpaware_ctr/seqs/{celltype}_nfrqtl_open.fasta",
        close_merge="../../results/results_03062025/motif_meme_snpaware_ctr/seqs/{celltype}_nfrqtl_close.fasta"
    envmodules:
        "Bioinformatics",
        "samtools/1.13-fwwss5n"
    shell:
        """
        cat {params.open_dir}chr*_open.fasta > {output.open_merge}
        cat {params.close_dir}chr*_close.fasta > {output.close_merge}
        samtools faidx {output.open_merge}
        samtools faidx {output.close_merge}
        """

rule sea:
    input:
        qtl="../../results/results_10102024/motif_meme_qtl/seqs/{celltype}_{region}qtl.fasta"
    output:
        "../../results/results_12032024/motif_meme_qtl/sea_{celltype}_{region}/sea.tsv"
    envmodules:
        "Bioinformatics",
        "meme/5.5.5"
    params:
        outdir="../../results/results_12032024/motif_meme_qtl/sea_{celltype}_{region}/",
        motif="../../data/motif/H12CORE_meme_format.meme",
        celltype="{celltype}",
        region="{region}"

    shell:
        """
        sea --p {input.qtl} --m {params.motif} --order 2 --thresh 1443 --seed 12 --oc {params.outdir}
        """

rule sea_coloc:
    input:
        nfrcoloc="../../results/results_02182025/motif_meme_coloc/seqs/{celltype}_nfrqtl_coloc.fasta",
        nfrno="../../results/results_02182025/motif_meme_coloc/seqs/{celltype}_nfrqtl_no_coloc.fasta",
    output:
        nfrcoloc="../../results/results_02182025/motif_meme_coloc/sea_{celltype}/sea.tsv"
        #nfrno="../../results/results_10102024/motif_meme_coloc/sea_{celltype}_nfrno/sea.tsv"
    envmodules:
        "Bioinformatics",
        "meme/5.5.5"
    params:
        outdir_nfrcoloc="../../results/results_02182025/motif_meme_coloc/sea_{celltype}/",
        #outdir_nfrno="../../results/results_10102024/motif_meme_coloc/sea_{celltype}_nfrno/",
        motif="../../data/motif/H12CORE_meme_format.meme",
    shell:
        """
        sea --p {input.nfrcoloc} --m {params.motif} --n {input.nfrno} --thresh 1443 --seed 12 --oc {params.outdir_nfrcoloc}

        """

rule sea_cit:
    input:
        "../../results/results_04232025/motif_meme_cit/seqs/{celltype}_cit.fasta"
    output:
         "../../results/results_04232025/motif_meme_cit/sea_{celltype}/sea.tsv"
    envmodules:
        "Bioinformatics",
        "meme/5.5.5"
    params:
        outdir="../../results/results_04232025/motif_meme_cit/sea_{celltype}/s",
        #outdir_nfrno="../../results/results_10102024/motif_meme_coloc/sea_{celltype}_nfrno/",
        motif="../../data/motif/H12CORE_meme_format.meme",
    shell:
        """
        sea --p {input} --m {params.motif} --order 2 --thresh 1443 --seed 12 --oc {params.outdir}

        """

rule sea_cit_merged:
    input:
        "../../results/results_04232025/motif_meme_cit/seqs/Type_1_cit.fasta"
    output:
         "../../results/results_04232025/motif_meme_cit/sea_merged/sea.tsv"
    envmodules:
        "Bioinformatics",
        "meme/5.5.5"
    params:
        indir="../../results/results_04232025/motif_meme_cit/seqs/",
        outdir="../../results/results_04232025/motif_meme_cit/sea_merged/",
        #outdir_nfrno="../../results/results_10102024/motif_meme_coloc/sea_{celltype}_nfrno/",
        motif="../../data/motif/H12CORE_meme_format.meme",
    shell:
        """
        cat {params.indir}*_cit.fasta > {params.outdir}combined.fasta
        sea --p {params.outdir}combined.fasta --m {params.motif} --order 2 --thresh 1443 --seed 12 --oc {params.outdir}

        """

rule sea_peaks:
    input:
        peak="../../results/results_11252024/motif_meme_peak/{celltype}_{region}Peak.fasta"
    output:
        "../../results/results_12032024/motif_meme_peak/sea_{celltype}_{region}Peak/sea.tsv"
    envmodules:
        "Bioinformatics",
        "meme/5.5.5"
    params:
        outdir="../../results/results_12032024/motif_meme_peak/sea_{celltype}_{region}Peak/",
        motif="../../data/motif/H12CORE_meme_format.meme",
        celltype="{celltype}",
        region="{region}"

    shell:
        """
        sea --p {input.peak} --m {params.motif} --order 2 --thresh 1443 --seed 12 --oc {params.outdir}
        """

rule sea_snpaware:
    input:
        open_merge="../../results/results_01032025/motif_meme_snpaware/seqs/{celltype}_nfrqtl_open.fasta",
        close_merge="../../results/results_01032025/motif_meme_snpaware/seqs/{celltype}_nfrqtl_close.fasta"
    output:
        "../../results/results_01032025/motif_meme_snpaware/sea_{celltype}/sea.tsv"
    envmodules:
        "Bioinformatics",
        "meme/5.5.5"
    params:
        outdir="../../results/results_01032025/motif_meme_snpaware/sea_{celltype}/",
        motif="../../data/motif/H12CORE_meme_format.meme"
    shell:
        """
        sea --p {input.open_merge} --m {params.motif} --n {input.close_merge} --pvalue 0.05 --seed 12 --oc {params.outdir}
        """

rule sea_snpaware_ctr:
    input:
        open_merge="../../results/results_03062025/motif_meme_snpaware_ctr/seqs/{celltype}_nfrqtl_open.fasta",
        close_merge="../../results/results_03062025/motif_meme_snpaware_ctr/seqs/{celltype}_nfrqtl_close.fasta"
    output:
        "../../results/results_03062025/motif_meme_snpaware_ctr/sea_{celltype}/sea.tsv"
    envmodules:
        "Bioinformatics",
        "meme/5.5.5"
    params:
        outdir="../../results/results_03062025/motif_meme_snpaware_ctr/sea_{celltype}/",
        motif="../../data/motif/H12CORE_meme_format.meme"
    shell:
        """
        sea --p {input.open_merge} --m {params.motif} --n {input.close_merge} --pvalue 0.05 --seed 12 --oc {params.outdir}
        """