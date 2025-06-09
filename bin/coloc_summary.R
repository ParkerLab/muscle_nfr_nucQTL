#!/usr/bin/env Rscript

# libraries
options(stringsAsFactors=FALSE)
library(coloc)
library(dplyr)
library(stringr)
library(susieR)
library(optparse)

option_list <- list(
    make_option(c("--rda_dir"), type = "character", help = "[Required] Susie model and coloc directory"),
    make_option(c("--GWAS_annot"), type = "character", help = "[Required] Checking which GWAS summary stats?"),
    make_option(c("--filter_window"), type = "numeric", help = "[Required] the mid-point distance between the peaks"),
    make_option(c("--outdir"), type = "character", help = "[Required] coloc summary results output directory")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

df = as.data.frame(matrix(nrow=0,ncol=13))
for (celltype in c('Type_1', 'Type_2a', 'Type_2x', 'Endothelial', 'Mesenchymal_Stem_Cell')){
    print(celltype)
    summary_path = paste0(opts$rda_dir, celltype, '_nfr', '/')
    print(summary_path)
    print(nrow(list.files(path=summary_path, pattern='*summary_100000win.Rda')))
    for (file in list.files(path=summary_path, pattern='*summary_100000win.Rda')){
        load(paste0(summary_path, file))
        if (!is.null(df_res_nfr) && nrow(df_res_nfr) > 0){
            df_res_nfr$Celltype=celltype
            #print(file)
            df = rbind(df, df_res_nfr)
        }
        #if (nrow(df_res) > 0){
        #    df_res$Celltype = celltype
        #    df = rbind(df, df_res)
        #}
    }
}

write.table(df, paste0(opts$outdir, 'nfrqtls_', opts$GWAS_annot, '_gwas_coloc.tsv'), quote=F)