#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)
library(locuscomparer)
library(dplyr)
library(susieR)
library(tidyr)
library(ggplot2)
library(gridExtra)

option_list <- list(
    make_option(c("--celltype"), type = "character", help = "[Required] QTLs from which celltype?"),
    make_option(c("--outdir_coloc"), type = "character", help = "[Required] Output directory for coloc results"),
    make_option(c("--outdir_plot"), type = "character", help = "[Required] Where to save the output plots?"),
    make_option(c("--outdir_susie"), type = "character", help = "[Required] Directory of the susie models"),
    make_option(c("--filter_window"), type = "numeric", help = "[Required] filter_window size for coloc (bp)"),
    make_option(c("--PPH4_threshold"), type = "numeric", help = "[Required] Filtering threshold for PPh4 value"),
    make_option(c("--nfr_count"), type = "character", help = "[Required] Directory of the normalized nfr count matrix"),
    make_option(c("--nuc_count"), type = "character", help = "[Required] Directory of the normalized nuc count matrix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

celltype=opts$celltype
outdir_coloc=opts$outdir_coloc
outdir_plot=opts$outdir_plot
outdir_susie=opts$outdir_susie
filter_window=opts$filter_window
PPH4_threshold=opts$PPH4_threshold
#-------------------plot functions-------------------
plot_scatter <- function(nfr_count_mat, nuc_count_mat, nfr_peak_id, nuc_peak_id, pph4, outfilename){
    # get data points
    nfr_col = nfr_count_mat %>% as.data.frame() %>% select(all_of(nfr_peak_id))
    nuc_col = nuc_count_mat %>% as.data.frame() %>% select(all_of(nuc_peak_id))
    # get sample id
    nfr_names = nfr_col %>% rownames()
    nfr_names = sapply(strsplit(nfr_names, "__"), function(x) x[[1]])
    nuc_names = nuc_col %>% rownames()
    nuc_names = sapply(strsplit(nuc_names, "__"), function(x) x[[1]])
    # add names col
    rownames(nfr_col) = nfr_names
    rownames(nuc_col) = nuc_names
    plot_df = merge(nfr_col, nuc_col, by = 'row.names', all = TRUE)
    # rename cols
    colnames(plot_df) = c("Sample", "NFR_peak", "Nuc_peak")
    # plot                   
    g = plot_df %>% ggplot(aes(x=NFR_peak, y=Nuc_peak)) +
        geom_point() +
        geom_smooth(method=lm) +
        ggtitle(paste("Nucleosomal peak vs. NFR peak accessibility", "PPH4:", round(pph4, 4))) +
        xlab(paste("NFR peak accessibility (normalized count)", nfr_peak_id)) + ylab(paste("Nucleosomal peak accessibility (normalized count)", nuc_peak_id))
    # save the plot
    ggsave(outfilename)
}
plot_locuscompare <- function(nfr_df, nuc_df, hit_nfr, hit_nuc, pph4, outfilename){
    # Determine the range for the x-axis based on the combined data
    x_min <- min(c(nfr_df$position, nuc_df$position))
    x_max <- max(c(nfr_df$position, nuc_df$position))
    # nfr locuscompare plot
    g_nfr = nfr_df %>% ggplot(aes(x=position, y=-log10(pval), label=rsid)) + 
    geom_point(aes(colour = r2), size=2) +
    ggtitle(paste("NFR lead snp:", hit_nfr, "PPH4:", round(pph4, 4))) +
    xlab(element_blank()) +
    ylab("nfrQTL - log10(P)") +
    scale_colour_steps2() +
    theme_classic() +
    guides(color = "none", shape = "none") +
    geom_text(data=subset(nfr_df, rsid == hit_nfr), nudge_x=1000, nudge_y=0.5) +
    geom_point(aes(x=nfr_df[nfr_df$rsid == hit_nfr,]$position,y=-log10(nfr_df[nfr_df$rsid == hit_nfr,]$pval)),
               colour="red", size=2) + xlim(x_min, x_max)
    # nuc locuscompare plot
    g_nuc = nuc_df %>% ggplot(aes(x=position, y=-log10(pval), label=rsid)) + 
    geom_point(aes(colour = r2), size=2) +
    ggtitle(paste("Nucleosomal lead snp:", hit_nuc)) +
    xlab(paste(chr, "(bp)")) + 
    ylab("nucQTL - log10(P)") +
    scale_colour_steps2() +
    theme_classic() +
    guides(color = "none", shape = "none") +
    geom_text(data=subset(nuc_df, rsid == hit_nuc), nudge_x=1000, nudge_y=0.5) +
    geom_point(aes(x=nuc_df[nuc_df$rsid == hit_nuc,]$position,y=-log10(nuc_df[nuc_df$rsid == hit_nuc,]$pval)),
               colour="red", size=2) + xlim(x_min, x_max)
    # compare nfr with nuc as a scatter plot
    merge_df = merge(nfr_df, nuc_df, by=c('rsid', 'position'), suffixes=c("_nfr", "_nuc"))
    g_cor = merge_df %>% ggplot(aes(x=-log10(pval_nfr), y=-log10(pval_nuc), label=rsid)) + 
    geom_point(aes(colour = r2_nfr), size=2) +
    xlab("nfrQTL - log10(P)") + 
    ylab("nucQTL - log10(P)") +
    scale_colour_steps2() +
    theme_classic() +
    geom_text(data=subset(merge_df, rsid == hit_nfr), nudge_x=0.5, nudge_y=0.5) +
    geom_point(aes(x=-log10(merge_df[merge_df$rsid == hit_nfr,]$pval_nfr), y=-log10(merge_df[merge_df$rsid == hit_nfr,]$pval_nuc)), colour="red", size=2) +
    guides(x=guide_axis(n.dodge=2))
    # arrange and save
    g = arrangeGrob(g_cor, g_nfr, g_nuc, layout_matrix=matrix(c(2,2,1,1,3,3,1,1), byrow = TRUE, nrow=2))
    ggsave(file=outfilename, g, width = 10, height = 5, dpi = 300)
}
#-------------------process coloc df-------------------
df_res_all = as.data.frame(matrix(nrow=0,ncol=12))
for (i in 1:22){
    chr = paste0('chr', i)
    # load the summary
    load(paste0(outdir_coloc, "coloc_", celltype, "/", chr, "/", chr, "_res_summary_", filter_window, "win.Rda"))
    df_res_all = rbind(df_res_all, df_res)
}
# trait df - now more like celltype df
trait_df = df_res_all
# get peak id col NFRPeak_NucleoPeak
trait_df$peak_id = paste(trait_df$NFR_peak, trait_df$Nucleosomal_peak, sep="_")
# add a chr col
trait_df$Chr = apply(trait_df, 1, FUN=function(x) strsplit(x['NFR_peak'], split="[.]")[[1]][1])
# save all results
output_filename = paste0(outdir_plot, "coloc_plots_", celltype, "/coloc_summary_", celltype, "_", filter_window, ".tsv")
write.table(trait_df, 
            file=output_filename, 
            sep='\t', quote=F)
#-------------------prepare for plotting-------------------
# filter and plot
trait_df = trait_df %>% filter(PP.H4.abf > PPH4_threshold)
## count mat locations
# HMMRATAC_nfr_mat = "../../data/normalized_counts/HMMRATAC_NFR_counts_multiotsu_normalized.rds"
# HMMRATAC_nuc_mat = "../../data/normalized_counts/HMMRATAC_Nucleosomal_counts_multiotsu_normalized.rds"
# NucleoATAC_nfr_mat = "../../data/normalized_counts/NuclaoATAC_NFR_counts_multiotsu_byreads_normalized.rds"
# NucleoATAC_nuc_mat = "../../data/normalized_counts/NucleoATAC_Nucleosomal_counts_multiotsu_200ext_04282023_filtered_normalized.rds"
# if (trait == "HMMRATAC"){
#         nfr_count_mat = readRDS(HMMRATAC_nfr_mat)
#         nuc_count_mat = readRDS(HMMRATAC_nuc_mat)
# } else if (trait == "NucleoATAC"){
#         nfr_count_mat = readRDS(NucleoATAC_nfr_mat)
#         nuc_count_mat = readRDS(NucleoATAC_nuc_mat)
# }
# read in the count mat
nfr_count_mat = readRDS(opts$nfr_count)
nuc_count_mat = readRDS(opts$nuc_count)
                     
# plot here
for (i in 1:nrow(trait_df)){
    #-------------------locuscompare-------------------
    # specify chr
    chr = trait_df[i,]$Chr
    file_nfr = paste0(outdir_susie, celltype, "_nfr_susie/", chr, "/",trait_df[i,]$NFR_peak, "_susie_dat.Rda")
    file_nuc = paste0(outdir_susie, celltype, "_nuc_susie/", chr, "/",trait_df[i,]$Nucleosomal_peak, "_susie_dat.Rda")
    # lead snp info
    hit_nfr = trait_df[i,]$hit1
    hit_nuc = trait_df[i,]$hit2
    # pph4 info
    pph4 = trait_df[i,]$PP.H4.abf
    # load nfr coloc dat
    load(file=file_nfr)
    nfr_dat = coloc_d1
    LD_nfr = nfr_dat$LD %>% as.data.frame() %>% select(all_of(hit_nfr)) %>% c()
    # square it?
    r2_nfr = nfr_dat$LD %>% as.data.frame() %>% select(all_of(hit_nfr)) %>% '^'(2) %>% c()
    nfr_df = data.frame(nfr_dat$snp, nfr_dat$pvalues, nfr_dat$position, LD_nfr, r2_nfr)
    names(nfr_df) = c('rsid', 'pval', 'position', 'LD', 'r2')
    nfr_df$lead = ifelse(nfr_df$rsid == hit_nfr, 1, 2)
    # load nuc coloc_dat
    load(file=file_nuc)
    nuc_dat = coloc_d1
    LD_nuc = nuc_dat$LD %>% as.data.frame() %>% select(all_of(hit_nuc)) %>% c()
    r2_nuc = nuc_dat$LD %>% as.data.frame() %>% select(all_of(hit_nuc)) %>% '^'(2) %>% c()
    nuc_df = data.frame(nuc_dat$snp, nuc_dat$pvalues, nuc_dat$position, LD_nuc, r2_nuc)
    names(nuc_df) = c('rsid', 'pval', 'position', 'LD', 'r2')
    nuc_df$lead = ifelse(nuc_df$rsid == hit_nuc, 1, 2)
    #-------------------locuscompare plot------------------
    peak_id = trait_df[i,]$peak_id
    outfilename_locuscompare = paste0(outdir_plot, "coloc_plots_", celltype, "/locuscomapreplot/locuscomapreplot_", peak_id, "_", celltype, "_",  filter_window, "win_PPH4", PPH4_threshold, "_", i, ".png")
    plot_locuscompare(nfr_df, nuc_df, hit_nfr, hit_nuc, pph4, outfilename_locuscompare)
    #-------------------scatter plot------------------
    chr = trait_df[i,]$Chr
    nfr_peak_id = trait_df[i,]$NFR_peak
    nuc_peak_id = trait_df[i,]$Nucleosomal_peak
    outfilename_scatter = paste0(outdir_plot, "coloc_plots_", celltype, "/scatterplot/scatterplot_", peak_id, "_", celltype, "_",  filter_window, "win_PPH4", PPH4_threshold, "_", i, ".png")
    plot_scatter(nfr_count_mat, nuc_count_mat, nfr_peak_id, nuc_peak_id, pph4, outfilename_scatter)
}
#trait_df$peak_id = paste(trait, trait_df$NFR_peak, trait_df$Nucleosomal_peak, sep="_")
# save output
output_filename = paste0(outdir_plot, "coloc_plots_", celltype, "/coloc_summary_", celltype, "_", filter_window, "win_PPH4", PPH4_threshold, ".tsv")
write.table(trait_df, 
            file=output_filename, 
            sep='\t', quote=F)

