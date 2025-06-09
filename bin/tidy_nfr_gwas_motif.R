#!/usr/bin/env Rscript

# libraries
options(stringsAsFactors=FALSE)
library(dplyr)
library(optparse)

option_list <- list(
    make_option(c("--nfrqtl"), type = "character", help = "[Required] nfrqtl permutation pass result"),
    make_option(c("--gwas"), type = "character", help = "[Required] nfr-gwas colocalization results"),
    make_option(c("--celltype"), type = "character", help = "[Required] which celltype"),
    make_option(c("--outdir"), type = "character", help = "[Required] outdir to save motif peaks"),
    make_option(c("--trait"), type = "character", help = "[Required] which gwas trait")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)
# celltype
celltype = opts$celltype
# read in the qtl results
nfrqtl = read.csv(opts$nfrqtl)
# read in the coloc results - qtl
gwas = read.csv(opts$gwas, sep = '\t')

# filter for significant qtls
sig_nfrqtl = nfrqtl %>% filter(adj_beta_qval < 0.05)

# filter for coloc qtls
gwas_sig = gwas %>% filter(PP.H4.abf > 0.5) %>% filter(Celltype == celltype) %>% select(c('hit1', 'PP.H4.abf', 'QTL_peak', 'gwas_peak'))

# motif peaks
# motif peaks
sig_nfrqtl %>% 
    merge(gwas_sig, by.x=c('phe_id'), by.y=c('QTL_peak'), all.x = TRUE) %>% 
    filter((!is.na(gwas_peak))) %>%
    select(c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd')) %>%
    write.table(file=paste0(opts$outdir, celltype, '_nfr_coloc_gwas_', opts$trait, '.homerpeak'), sep = '\t', quote=F, row.names = F, col.names = F)

# background peaks
sig_nfrqtl %>% 
    merge(gwas_sig, by.x=c('phe_id'), by.y=c('QTL_peak'), all.x = TRUE) %>% 
    filter((is.na(gwas_peak))) %>%
    select(c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd')) %>%
    write.table(file=paste0(opts$outdir, celltype, '_nfr_coloc_gwas_', opts$trait, '_bkg.homerpeak'), sep = '\t', quote=F, row.names = F, col.names = F)