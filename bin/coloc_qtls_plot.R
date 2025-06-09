#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(dplyr)
library(locuscomparer)
library(susieR)
library(tidyr)
library(optparse)

option_list <- list(
    make_option(c("--trait"), type = "character", help = "[Required] Which qtl method comparing?"),
    make_option(c("--outdir_tidy"), type = "character", help = "[Required] Directory of the coloc tidy results"),
    make_option(c("--outdir_plot"), type = "character", help = "[Required] Where to save the output plots?"),
    make_option(c("--outdir_susie"), type = "character", help = "[Required] Directory of the susie models"),
    make_option(c("--filter_window"), type = "numeric", help = "[Required] filter_window size for coloc (bp)"),
    make_option(c("--PPH4_threshold"), type = "numeric", help = "[Required] Filtering threshold for PPh4 value"),
    make_option(c("--chr"), type = "character", help = "[Required] Which chr plotting?"),
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

trait=opts$trait
outdir_tidy=opts$outdir_tidy
outdir_plot=opts$outdir_plot
outdir_susie=opts$outdir_susie
filter_window=opts$filter_window
PPH4_threshold=opts$PPH4_threshold
chr=opts$chr

# process input
trait_df = read.table(paste0(outdir_tidy, "coloc_plots_", trait, "/coloc_summary_HMMRATAC_",filter_window, "win_PPH4",PPH4_threshold, ".tsv"))
trait_df$Chr = apply(trait_df, 1, FUN=function(x) strsplit(x['NFR_peak'], split="[.]")[[1]][1])
# filter to chr
trait_df_chr = trait_df %>% filter(Chr == chr)     

for (i in 1:nrow(trait_df_chr)){
    file_nfr = paste0(outdir_susie, trait, "_NFR_sigQTL_susie_models/", chr, "/",trait_df_chr[i,]$NFR_peak, "_susie_dat.Rda")
    file_nuc = paste0(outdir_susie, trait, "_Nucleosomal_sigQTL_susie_models/", chr, "/",trait_df_chr[i,]$Nucleosomal_peak, "_susie_dat.Rda")
    # load nfr coloc dat
    load(file=file_nfr)
    nfr_df = data.frame(coloc_d1$snp, coloc_d1$pvalues)
    names(nfr_df) = c('rsid', 'pval')
    # load nuc coloc_dat
    load(file=file_nuc)
    nuc_df = data.frame(coloc_d1$snp, coloc_d1$pvalues)
    names(nuc_df) = c('rsid', 'pval')
    print(file_nfr)
    print(file_nuc)
    # locuscompare plot
    locuscompare_filename = paste0(outdir_plot, chr, '/', trait_df_chr[i,]$peak_id, '_locuscompare.png')
    png(locuscompare_filename)
    locuscompare(in_fn1 = nfr_df,
             in_fn2 = nuc_df,
             title = 'nfrQTL', title2 = 'nucQTL', genome="hg38")
    dev.off()
}