#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)
library(dplyr)
library(tidyverse)

option_list <- list(
    make_option(c("--celltype"), type = "character", help = "[Required] Which celltype are we looking at"),
    make_option(c("--region"), type = "character", help = "[Required] Which region? nfr or nuc?"),
    make_option(c("--input_dir"), type = "character", help = "[Required] Directory that stored the bw scores?"),
    make_option(c("--output"), type = "character", help = "[Required] Where to store the tidied count matrix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# input
celltype = opts$celltype
region = opts$region
input_dir = opts$input_dir
output = opts$output

# total number of samples - should be 284
n = length(list.files(input_dir))
init_df = read.table(paste0(input_dir, list.files(input_dir)[1]), 
                     sep = '\t', 
                     col.names=c('name', 'size', 'covered', 'sum', 'mean0', 'mean', 'min', 'max'))
init_df = init_df %>% select(name) %>% filter(!str_detect(name, 'chrX|chrY'))
# set up the two quantification df
df_max = init_df
df_mean = init_df
for(i in 1:n){
    #print(i)
    summary_path = paste0(input_dir, list.files(input_dir)[i])
    # colname
    colname = strsplit(list.files(input_dir)[i], '_')[[1]][1]
    # reformat colname
    colname = paste0(gsub('-', '_', colname), '_', region, '_bw')
    # read in the df
    df = read.table(summary_path, 
                    sep = '\t', 
                    col.names=c('name', 'size', 'covered', 'sum', 'mean0', 'mean', 'min', 'max'))
    # extract scores
    df_max_tmp = df %>% select(name, max)
    df_mean_tmp = df %>% select(name, mean)
    # merge
    df_max = merge(df_max, df_max_tmp)
    df_mean = merge(df_mean, df_mean_tmp)
    # rename
    names(df_max)[names(df_max) == 'max'] <- colname
    names(df_mean)[names(df_mean) == 'mean'] <- colname
}

# set index
rownames(df_max) <- df_max$name
df_max$name <- NULL

rownames(df_mean) <- df_mean$name
df_mean$name <- NULL

#convert to matrix
m_max = as.matrix(df_max)
m_mean = as.matrix(df_mean)
#save as rds
saveRDS(m_max, file=paste0(output, 'max_tidy.rds'))
saveRDS(m_mean, file=paste0(output, 'mean_tidy.rds'))