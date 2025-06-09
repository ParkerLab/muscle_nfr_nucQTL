#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)
library(coloc)
library(dplyr)
library(stringr)
library(susieR)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(patchwork)

option_list <- list(
    make_option(c("--LD_file"), type = "character", help = "[Required] Which gwas lead var is examining?"),
    make_option(c("--celltype"), type = "character", help = "[Required] Which celltype?"),
    make_option(c("--outdir_coloc"), type = "character", help = "[Required] Where is the directory for coloc results"),
    make_option(c("--GWAS_annot"), type = "character", help = "[Required] Which GWAS summary stats?"),
    make_option(c("--region"), type = "character", help = "[Required] which type of QTL? NFR vs. Nucleosoaml"),
    make_option(c("--PPH4_threshold"), type = "numeric", help = "[Required] Filtering threshold for PPh4 value"),
    make_option(c("--outdir_plot"), type = "character", help = "[Required] Where to save the output plots?"),
    make_option(c("--dir_qtl_susie"), type = "character", help = "[Required] Directory of the QTL susie models"),
    make_option(c("--dir_gwas_susie"), type = "character", help = "[Required] Directory of the GWAS susie models"),
    make_option(c("--filter_window"), type = "numeric", help = "[Required] filter_window size for coloc (bp)"),
    make_option(c("--gwas_summary"), type = "character", help = "[Required] Summary file per gwas signal coloc results")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)
#-------------------plot functions-------------------
plot_locuscompare <- function(qtl_df, gwas_df, hit_qtl, hit_gwas, region, pph4, outfilename){
    # Determine the range for the x-axis based on the combined data
    x_min <- min(c(qtl_df$position, gwas_df$position))
    x_max <- max(c(qtl_df$position, gwas_df$position))
    # qtl locuscompare plot
    g_qtl = qtl_df %>% ggplot(aes(x=position, y=-log10(pval), label=rsid)) + 
    geom_point(aes(colour = r2), size=2) +
    ggtitle(paste("QTL lead snp -", region, hit_qtl, "PPH4", round(pph4, 4))) +
    xlab(element_blank()) +
    ylab(paste0(region, "QTL - log10(P)")) +
    scale_colour_steps2() +
    theme_classic() +
    guides(color = "none", shape = "none") +
    geom_text(data=subset(qtl_df, rsid == hit_qtl), nudge_x=1000, nudge_y=0.5) +
    geom_point(aes(x=qtl_df[qtl_df$rsid == hit_qtl,]$position,y=-log10(qtl_df[qtl_df$rsid == hit_qtl,]$pval)),
            colour="red", size=2) + xlim(x_min, x_max)
    # gwas locuscompare plot
    g_gwas = gwas_df %>% ggplot(aes(x=position, y=-log10(pval), label=rsid)) + 
    geom_point(aes(colour = r2), size=2) +
    ggtitle(paste("GWAS lead snp -", hit_gwas)) +
    xlab(element_blank()) +
    ylab("GWAS - log10(P)") +
    scale_colour_steps2() +
    theme_classic() +
    guides(color = "none", shape = "none") +
    geom_text(data=subset(gwas_df, rsid == hit_gwas), nudge_x=1000, nudge_y=0.5) +
    geom_point(aes(x=gwas_df[gwas_df$rsid == hit_gwas,]$position,y=-log10(gwas_df[gwas_df$rsid == hit_gwas,]$pval)),
            colour="red", size=2) + xlim(x_min, x_max)
    # compare qtl with gwas as a scatter plot
    merge_df = merge(qtl_df, gwas_df, by=c('rsid', 'position'), suffixes=c("_qtl", "_gwas"))
    g_cor = merge_df %>% ggplot(aes(x=-log10(pval_qtl), y=-log10(pval_gwas), label=rsid)) + 
    geom_point(aes(colour = r2_gwas), size=2) +
    xlab(paste(region, "QTL - log10(P)")) + 
    ylab("GWAS - log10(P)") +
    scale_colour_steps2() +
    theme_classic() +
    geom_text(data=subset(merge_df, rsid == hit_gwas), nudge_x=0.5, nudge_y=0.5) +
    geom_point(aes(x=-log10(merge_df[merge_df$rsid == hit_gwas,]$pval_qtl), y=-log10(merge_df[merge_df$rsid == hit_gwas,]$pval_gwas)), colour="red", size=2) +
    guides(x=guide_axis(n.dodge=2))
    # arrange and save
    g = arrangeGrob(g_cor, g_qtl, g_gwas, layout_matrix=matrix(c(2,2,1,1,3,3,1,1), byrow = TRUE, nrow=2))
    # arrange and save using patchwork
    #combined_plot <- (g_qtl / g_gwas) | g_cor
    ggsave(file=outfilename, g, width = 10, height = 5, dpi = 300)
}

# input variables
LD_file = opts$LD_file
# which celltype
celltype = opts$celltype
# coloc dir
outdir_coloc = opts$outdir_coloc
# GWAS annotation
GWAS_annot = opts$GWAS_annot
# region
region = opts$region
# filtering threshold
PPH4_threshold = opts$PPH4_threshold
# output plots dir
outdir_plot = opts$outdir_plot
# QTL susie dir
dir_qtl_susie = opts$dir_qtl_susie
# GWAS susie dir
dir_gwas_susie = opts$dir_gwas_susie
# GWAS dat input filepath
gwas_summary = opts$gwas_summary
# coloc window
filter_window = opts$filter_window
filter_window = format(filter_window, scientific = FALSE)

# load the summary file
load(gwas_summary)
# if (region == "NFR"){
#     df_res = df_res_nfr
# }
# if (region == "Nucleosomal"){
#     df_res = df_res_nuc
# }
df_res = df_res_nfr
# loop and plot
if (nrow(df_res) > 0){
    # assign chr
    df_res$Chr = apply(df_res, 1, FUN=function(x) strsplit(x['gwas_peak'], split="-")[[1]][1])
    # filter by PPH4
    df_res = df_res %>% filter(PP.H4.abf > PPH4_threshold)
    # if still greater than 0
    if (nrow(df_res) > 0){
        for (i in 1:nrow(df_res)){
            chr = df_res[i,]$Chr
            file_qtl = paste0(dir_qtl_susie, celltype, "_", region, "_susie/", chr, "/", df_res[i,11], "_susie_dat.Rda")
            file_gwas = paste0(dir_gwas_susie, GWAS_annot, "_susie_models/", LD_file, "_susie_dat.Rda")
            # lead snp info = wrong! reversed!
            hit_qtl = df_res[i,]$hit2
            hit_gwas = df_res[i,]$hit1
            # add PPH4
            pph4 = df_res[i,]$PP.H4.abf
            # load qtl coloc dat
            load(file=file_qtl)
            qtl_dat = coloc_d1
            LD_qtl = qtl_dat$LD %>% as.data.frame() %>% select(all_of(hit_qtl)) %>% c()
            # square it?
            r2_qtl = qtl_dat$LD %>% as.data.frame() %>% select(all_of(hit_qtl)) %>% '^'(2) %>% c()
            qtl_df = data.frame(qtl_dat$snp, qtl_dat$pvalues, qtl_dat$position, LD_qtl, r2_qtl)
            names(qtl_df) = c('rsid', 'pval', 'position', 'LD', 'r2')
            qtl_df$lead = ifelse(qtl_df$rsid == hit_qtl, 1, 2)
            # load gwas coloc_dat
            load(file=file_gwas)
            gwas_dat = coloc_d1
            LD_gwas = gwas_dat$LD %>% as.data.frame() %>% select(all_of(hit_gwas)) %>% c()
            r2_gwas = gwas_dat$LD %>% as.data.frame() %>% select(all_of(hit_gwas)) %>% '^'(2) %>% c()
            gwas_df = data.frame(gwas_dat$snp, gwas_dat$pvalues, gwas_dat$position, LD_gwas, r2_gwas)
            names(gwas_df) = c('rsid', 'pval', 'position', 'LD', 'r2')
            gwas_df$lead = ifelse(gwas_df$rsid == hit_gwas, 1, 2)
            print(file_qtl)
            print(file_gwas)
            # plot
            ## filename
            filename_locuscompare = paste0(outdir_plot, "gwas_coloc_plots_", celltype, "_", filter_window, "win", "/", GWAS_annot, "_", region,
                                "/locuscomapreplot/locuscomapreplot_", LD_file, "_", region, "_",
                                df_res[i,11], "win_PPH4", PPH4_threshold, "_", i, ".png")
            ## locuscompare plot
            #plot_locuscompare(qtl_df, gwas_df, hit_qtl, hit_gwas, region, pph4, filename_locuscompare)
        }
    }
}
# summary df
df <- data.frame(Region = LD_file, num_pairs = nrow(df_res))
# save the table
write.table(df, 
            file=paste0(outdir_plot, "gwas_coloc_plots_", celltype, "_", filter_window, "win", "/", GWAS_annot, "_", region,
            "/locuscomapreplot/logs/", LD_file, "_win_PPH4", PPH4_threshold, ".tsv"), 
            sep='\t', quote=F)