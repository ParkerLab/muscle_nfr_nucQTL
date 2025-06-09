library(locuscomparer)
library(dplyr)
library(susieR)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(coloc)
library(tidyr)
library(UpSetR)

# read in coloc results
coloc_df <- read.table('../../results/results_05142024/nfr_nucqlts_coloc_05142024.tsv') %>% filter(!if_any(everything(), is.na))

# init track lists
l_n_nfr_cs_snps <- c()
l_n_nuc_cs_snps <- c()
l_n_nfr_cs_snps_inpeak <- c()
l_n_nuc_cs_snps_inpeak <- c()

susiedir = "../../results/results_05132024/"
nominaldir = "../../results/results_05102024/QTL_opt_results/"
for( i in 1:nrow(coloc_df)){
#for( i in 5000:5001){
    # get peak info
    nfrpeak = coloc_df[i,]$NFR_peak
    nucpeak = coloc_df[i,]$Nucleosomal_peak
    celltype = coloc_df[i,]$Celltype
    # get chr info
    chr = strsplit(nfrpeak, '\\.')[[1]][1]
    # get susie dir info
    nfr_cs_dir = paste0(susiedir, celltype, '_nfr_susie/', chr, '/', nfrpeak, '.susie_cset95.tsv')
    nuc_cs_dir = paste0(susiedir, celltype, '_nuc_susie/', chr, '/', nucpeak, '.susie_cset95.tsv')
    # read in
    nfr_cs_df <- read.table(nfr_cs_dir)
    nuc_cs_df <- read.table(nuc_cs_dir)
    # get total cs snps number
    n_nfr_cs_snps = nrow(nfr_cs_df)
    n_nuc_cs_snps = nrow(nuc_cs_df)
    ## append to list
    l_n_nfr_cs_snps <- c(l_n_nfr_cs_snps, n_nfr_cs_snps)
    l_n_nuc_cs_snps <- c(l_n_nuc_cs_snps, n_nuc_cs_snps)
    # get nominal qtl info
    nfr_nominal_dir = paste0(nominaldir, celltype, '_nfr/nominal_', chr, '.txt')
    nuc_nominal_dir = paste0(nominaldir, celltype, '_nuc/nominal_', chr, '.txt')
    # read in
    nfr_nominal_df <- read.table(nfr_nominal_dir, col.names=c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var', 'var_id', 
                                      'var_chr', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'slope', 'slope_se', 'best_hit'))
    nuc_nominal_df <- read.table(nuc_nominal_dir, col.names=c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var', 'var_id', 
                                      'var_chr', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'slope', 'slope_se', 'best_hit'))
    # filter by pheid
    nfr_peak_nominal_df <- nfr_nominal_df %>% filter(phe_id == nfrpeak)
    nuc_peak_nominal_df <- nuc_nominal_df %>% filter(phe_id == nucpeak)
    # merge with cs df
    nfr_merge <- merge(nfr_cs_df, nfr_peak_nominal_df, by.x='snp', by.y='var_id', how='left')
    nuc_merge <- merge(nuc_cs_df, nuc_peak_nominal_df, by.x='snp', by.y='var_id', how='left')
    # get numbers
    n_nfr_snps_inpeak <- nfr_merge %>% filter(dist_phe_var == 0) %>% nrow() %>% c()
    n_nuc_snps_inpeak <- nuc_merge %>% filter(dist_phe_var == 0) %>% nrow() %>% c()
    ## append to list
    l_n_nfr_cs_snps_inpeak <- c(l_n_nfr_cs_snps_inpeak, n_nfr_snps_inpeak)
    l_n_nuc_cs_snps_inpeak <- c(l_n_nuc_cs_snps_inpeak, n_nuc_snps_inpeak)
    # susie directory
    #coloc_df[i,]
}
coloc_df$n_nfr_cs_snps <- l_n_nfr_cs_snps
coloc_df$n_nuc_cs_snps <- l_n_nuc_cs_snps
coloc_df$n_nfr_cs_snps_inpeak <- l_n_nfr_cs_snps_inpeak
coloc_df$n_nuc_cs_snps_inpeak <- l_n_nuc_cs_snps_inpeak

write.table(coloc_df, '../../results/results_05142024/coloc_snps_inpeak.tsv', sep = '\t', quote=F)
