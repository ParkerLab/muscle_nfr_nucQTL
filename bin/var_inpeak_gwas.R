library(dplyr)
library(tidyr)
library(ggplot2)
library(Rfast)
library(optparse)

option_list <- list(
    make_option(c("--susiedir"), type = "character", help = "[Required] SuSiE results dir"),
    make_option(c("--chr"), type = "character", help = "[Required] Checking on which chromosome"),
    make_option(c("--celltype"), type = "character", help = "[Required] Checking on which chromosome"),
    make_option(c("--nfr_nominaldir"), type = "character", help = "[Required] nfrQTL nominal pass results"),
    make_option(c("--gwas_summary"), type = "character", help = "[Required] path for gwas-nfrQTL coloc results"),
    make_option(c("--output"), type = "character", help = "[Required] Var snps peak output directory")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# process input
chr = opts$chr
celltype = opts$celltype
susiedir = opts$susiedir
nfr_nominaldir = opts$nfr_nominaldir
# read in coloc df
coloc_df <- read.table(opts$gwas_summary) %>% filter(!if_any(everything(), is.na))
# filter to chr and celltype
coloc_df_tmp <- coloc_df %>% 
                    separate(QTL_peak, c("Chr", NA, NA), remove=FALSE) %>% 
                    filter(Chr == chr) %>%
                    filter(Celltype == celltype) %>%
                    filter(PP.H4.abf > 0.5)

# get nominal qtl info
nfr_nominal_dir = paste0(nfr_nominaldir, celltype, '_nfr/nominal_', chr, '.txt')
# read in
nfr_nominal_df <- read.table(nfr_nominal_dir, col.names=c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var', 'var_id', 
                                    'var_chr', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'slope', 'slope_se', 'best_hit'))

if(nrow(coloc_df_tmp) > 0){
    coloc_df_tmp$in_peak = 0
    for( i in 1:nrow(coloc_df_tmp)){
    #for( i in 5000:5001){
        # get peak info
        nfrpeak = coloc_df_tmp[i,]$QTL_peak
        gwaspeak = coloc_df_tmp[i,]$gwas_peak
        #celltype = coloc_df_tmp[i,]$Celltype
        # get chr info
        #chr = strsplit(nfrpeak, '\\.')[[1]][1]
        # get susie dir info
        nfr_cs_dir = paste0(susiedir, celltype, '_nfr_susie/', chr, '/', nfrpeak, '.susie_cset95.tsv')
        # read in
        nfr_cs_df <- read.table(nfr_cs_dir)
        # filter by pheid
        nfr_peak_nominal_df <- nfr_nominal_df %>% filter(phe_id == nfrpeak)
        # merge with cs df
        nfr_merge <- merge(nfr_cs_df, nfr_peak_nominal_df, by.x='snp', by.y='var_id', how='left')
        # filter for inpeak snps
        nfr_merge_tmp <- nfr_merge %>% filter(dist_phe_var == 0)
        # add inpeak col
        if (nfr_merge_tmp %>% nrow() > 0){
            coloc_df_tmp[i,]$in_peak = 1
        }
    }
}


write.table(coloc_df_tmp, opts$output, sep = '\t', quote=F)
