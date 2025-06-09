library(dplyr)
library(tidyr)
library(ggplot2)
library(Rfast)
library(optparse)

option_list <- list(
    make_option(c("--coloc_dir"), type = "character", help = "[Required] coloc results"),
    make_option(c("--susiedir"), type = "character", help = "[Required] SuSiE results dir"),
    make_option(c("--chr"), type = "character", help = "[Required] Checking on which chromosome"),
    make_option(c("--celltype"), type = "character", help = "[Required] Checking on which chromosome"),
    make_option(c("--nfr_nominaldir"), type = "character", help = "[Required] nfrQTL nominal pass results"),
    make_option(c("--nuc_nominaldir"), type = "character", help = "[Required] nucQTL nominal pass results"),
    make_option(c("--outdir"), type = "character", help = "[Required] Var snps peak output directory")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# process input
chr = opts$chr
celltype = opts$celltype
susiedir = opts$susiedir
nfr_nominaldir = opts$nfr_nominaldir
nuc_nominaldir = opts$nuc_nominaldir
# read in coloc df
coloc_df <- read.table(opts$coloc_dir) %>% filter(!if_any(everything(), is.na)) %>% filter(PP.H4.abf > 0.8)
# filter to chr and celltype
coloc_df_tmp <- coloc_df %>% 
                    separate(NFR_peak, c("Chr", NA, NA), remove=FALSE) %>% 
                    filter(Chr == chr) %>%
                    filter(Celltype == celltype)


# init track lists
l_n_nfr_cs_snps <- c()
l_n_nuc_cs_snps <- c()
l_n_nfr_cs_snps_inpeak <- c()
l_n_nuc_cs_snps_inpeak <- c()
row_list_nfr <- list()
row_list_nuc <- list()

# get nominal qtl info
nfr_nominal_dir = paste0(nfr_nominaldir, celltype, '_nfr/nominal_', chr, '.txt')
nuc_nominal_dir = paste0(nuc_nominaldir, celltype, '_nuc/nominal_', chr, '.txt')
# read in
nfr_nominal_df <- read.table(nfr_nominal_dir, col.names=c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var', 'var_id', 
                                    'var_chr', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'slope', 'slope_se', 'best_hit'))
nuc_nominal_df <- read.table(nuc_nominal_dir, col.names=c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var', 'var_id', 
                                    'var_chr', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'slope', 'slope_se', 'best_hit'))

for( i in 1:nrow(coloc_df_tmp)){
#for( i in 5000:5001){
    # get peak info
    nfrpeak = coloc_df_tmp[i,]$NFR_peak
    nucpeak = coloc_df_tmp[i,]$Nucleosomal_peak
    #celltype = coloc_df_tmp[i,]$Celltype
    # get chr info
    #chr = strsplit(nfrpeak, '\\.')[[1]][1]
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
    # filter by pheid
    nfr_peak_nominal_df <- nfr_nominal_df %>% filter(phe_id == nfrpeak)
    nuc_peak_nominal_df <- nuc_nominal_df %>% filter(phe_id == nucpeak)
    # merge with cs df
    nfr_merge <- merge(nfr_cs_df, nfr_peak_nominal_df, by.x='snp', by.y='var_id', how='left')
    nuc_merge <- merge(nuc_cs_df, nuc_peak_nominal_df, by.x='snp', by.y='var_id', how='left')
    # add inpeak col
    nfr_merge <- nfr_merge %>% mutate(in_peak = case_when(dist_phe_var == 0 ~ 1,
                             TRUE ~ 0))
    nuc_merge <- nuc_merge %>% mutate(in_peak = case_when(dist_phe_var == 0 ~ 1,
                             TRUE ~ 0))
    # get numbers
    n_nfr_snps_inpeak <- nfr_merge %>% filter(dist_phe_var == 0) %>% nrow() %>% c()
    n_nuc_snps_inpeak <- nuc_merge %>% filter(dist_phe_var == 0) %>% nrow() %>% c()
    ## append to list
    l_n_nfr_cs_snps_inpeak <- c(l_n_nfr_cs_snps_inpeak, n_nfr_snps_inpeak)
    l_n_nuc_cs_snps_inpeak <- c(l_n_nuc_cs_snps_inpeak, n_nuc_snps_inpeak)
    # susie directory
    #coloc_df_tmp[i,]
    # append df
    row_list_nfr[[i]] = nfr_merge
    row_list_nuc[[i]] = nuc_merge
}
coloc_df_tmp$n_nfr_cs_snps <- l_n_nfr_cs_snps
coloc_df_tmp$n_nuc_cs_snps <- l_n_nuc_cs_snps
coloc_df_tmp$n_nfr_cs_snps_inpeak <- l_n_nfr_cs_snps_inpeak
coloc_df_tmp$n_nuc_cs_snps_inpeak <- l_n_nuc_cs_snps_inpeak
# cs pip df
nfr_merge_chr <- do.call(rbind, row_list_nfr)
nuc_merge_chr <- do.call(rbind, row_list_nuc)

write.table(nfr_merge_chr, gsub('_info.tsv', '_nfr_cs_pip.tsv', opts$outdir), sep = '\t', quote=F)
write.table(nuc_merge_chr, gsub('_info.tsv', '_nuc_cs_pip.tsv', opts$outdir), sep = '\t', quote=F)
write.table(coloc_df_tmp, opts$outdir, sep = '\t', quote=F)
