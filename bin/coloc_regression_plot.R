#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(optparse)

option_list <- list(
    make_option(c("--coloc_results"), type = "character", help = "[Required] The path for qtl coloc results summary?"),
    make_option(c("--celltype"), type = "character", help = "[Required] which celltype?"),
    make_option(c("--output"), type = "character", help = "[Required] Output filename")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

coloc_df = read.table(opts$coloc_results)
celltype = opts$celltype

#nfr_qtl_path = paste0('../../results/results_05102024/QTL_opt_results/', celltype, '_nfr/merged_QTLresults_fdrcorr.csv')
#nuc_qtl_path = paste0('../../results/results_07162024/QTL_opt_results/', celltype, '_nuc/merged_QTLresults_fdrcorr.csv')
# select for fdr and select columns
#nfr_qtl_df = read.csv(nfr_qtl_path) %>% filter(adj_beta_qval < 0.05) %>% select(c('phe_id', 'slope'))
#nuc_qtl_df = read.csv(nuc_qtl_path) %>% filter(adj_beta_qval < 0.05) %>% select(c('phe_id', 'slope'))
# fitler for celltype
coloc_celltype = coloc_df %>% filter(Celltype == celltype)
# loop with nominal
coloc_slope_df <- data.frame()
for(i in 1:22){
    chr=paste0('chr', i)
    coloc_chr_df = coloc_celltype %>% filter(str_detect(NFR_peak, chr))
    print(chr)
    # read in nominal df
    nominal_nfr_df <- fread(paste0('../../results/results_05102024//QTL_opt_results//', celltype, '_nfr/nominal_', chr, '.txt'), sep = ' ',
                       colClasses = c('character', 'character', 'numeric', 'numeric', 
                                      'character', 'numeric', 'numeric', 'character',
                                      'character', 'numeric', 'numeric', 'numeric',
                                       'numeric', 'numeric', 'numeric', 'numeric'),
                       col.names=c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var', 'var_id', 
                                   'var_chr', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'slope', 'slope_se', 'best_hit'))
    nominal_nfr_df <- nominal_nfr_df %>% select(c('phe_id', 'var_id', 'slope', 'nom_pval'))
    #nuc
    nominal_nuc_df <- fread(paste0('../../results/results_07162024//QTL_opt_results//', celltype, '_nuc/nominal_', chr, '.txt'), sep = ' ',
                       colClasses = c('character', 'character', 'numeric', 'numeric', 
                                      'character', 'numeric', 'numeric', 'character',
                                      'character', 'numeric', 'numeric', 'numeric',
                                       'numeric', 'numeric', 'numeric', 'numeric'),
                       col.names=c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var', 'var_id', 
                                   'var_chr', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'slope', 'slope_se', 'best_hit'))
    nominal_nuc_df <- nominal_nuc_df %>% select(c('phe_id', 'var_id', 'slope', 'nom_pval'))
    # merge
    coloc_chr_df <- coloc_chr_df %>%
        filter(PP.H4.abf > 0.8) %>%
        merge(nominal_nfr_df, by.x=c('NFR_peak', 'hit1'), by.y=c('phe_id', 'var_id')) %>%
        rename('slope_nfr' = 'slope') %>%
        rename('nom_pval_nfr' = 'nom_pval') %>%
        merge(nominal_nuc_df, by.x=c('Nucleosomal_peak', 'hit1'), by.y=c('phe_id', 'var_id')) %>%
        rename('slope_nuc' = 'slope') %>%
        rename('nom_pval_nuc' = 'nom_pval')
    # append
    coloc_slope_df <- rbind(coloc_slope_df, coloc_chr_df)
}
# save results
write.table(coloc_slope_df, file=paste0(opts$output, '.tsv'), 
            sep='\t', quote=F)
coloc_slope_df %>% 
    ggplot(aes(x=slope_nuc, y=slope_nfr)) + 
    geom_point(color=4, alpha=0.1, size=3) +
    ggtitle(paste(celltype, 'regression plot (PPH4>0.8)'))
# save figure
ggsave(paste0(opts$output, '.png'))
