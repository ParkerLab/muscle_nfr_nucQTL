#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)
library(dplyr)
library(stringr)

option_list <- list(
    make_option(c("--celltype"), type = "character", help = "[Required] The celltype that we are checking"),
    make_option(c("--eqtl_dir"), type = "character", help = "[Required] eQTL results saved location"),
    make_option(c("--output"), type = "character", help = "[Required] eQTL summary saved output filename")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

eqtl_summary <- function(celltype, eqtl_dir){
    # file list
    eqtl_files = c()
    # list all folders
    eqtl_folders = list.files(eqtl_dir, pattern=celltype)
    for (folder in eqtl_folders){
        eqtl_folder_path = paste0(eqtl_dir, folder)
        # list all files
        eqtl_file_tmp = list.files(eqtl_folder_path, pattern="*.bed", full.names=T)
        # append
        eqtl_files = unlist(c(eqtl_files, eqtl_file_tmp))
    }
    # tracking
    # init tracking list
    gene_list = c()
    chr_list = c()
    pos_list = c()
    leadvar_list = c()
    path_list = c()
    # check all possible summarys stats
    for(file in eqtl_files){
        # process all information from filename
        eqtl_full_path = file
        file = strsplit(eqtl_full_path, split='/')[[1]][6]
        gene = strsplit(strsplit(strsplit(file, split=celltype)[[1]][2], '__')[[1]][1], split='\\.')[[1]][5]
        chr = strsplit(strsplit(strsplit(strsplit(file, split=celltype)[[1]][2], '__')[[1]][2], split='-')[[1]][2], split='_')[[1]][1]
        pos = strsplit(strsplit(strsplit(strsplit(file, split=celltype)[[1]][2], '__')[[1]][2], split='-')[[1]][2], split='_')[[1]][2]
        leadvar = strsplit(strsplit(strsplit(file, split=celltype)[[1]][2], '__')[[1]][2], split='-')[[1]][3] %>% str_replace('.cset.bed', '')
        # append information to list
        gene_list = c(gene_list, gene)
        chr_list = c(chr_list, chr)
        pos_list = c(pos_list, pos)
        leadvar_list = c(leadvar_list, leadvar)
        path_list = c(path_list, eqtl_full_path)
    }
    # conbime list
    my_list = list(gene = gene_list, 
                   chr = chr_list, 
                   pos = pos_list, 
                   leadvar = leadvar_list, 
                   path = path_list)
    # form df
    eqtl_summary_df = as.data.frame(my_list)
    
    return (eqtl_summary_df)
}

eqtl_summary_df = eqtl_summary(opts$celltype, opts$eqtl_dir)
write.table(eqtl_summary_df, 
            file=opts$output, 
            sep='\t', quote=F)