#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(qvalue)
library(stringr)

path = args[1]

df_NFR_x = read.table(path, sep = '\t', header=TRUE)
df_NFR_x$adj_beta_qval = qvalue(df_NFR_x$adj_beta_pval)$qvalues
df_NFR_x_significant = df_NFR_x[df_NFR_x$adj_beta_qval < 0.05,]
#paste0(path, ' Num of total peaks: ', nrow(df_NFR_x),'. Num of significant peaks: ', nrow(df_NFR_x_significant))
write.csv(df_NFR_x, str_replace(path,'.txt', '_fdrcorr.csv')) 
