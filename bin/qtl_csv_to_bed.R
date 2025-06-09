#!/usr/bin/env Rscript

# libraries
options(stringsAsFactors=FALSE)
library(dplyr)
library(optparse)

option_list <- list(
    make_option(c("--csv"), type = "character", help = "[Required] Input qtl results in .csv format (fdr-corrected)"),
    make_option(c("--out"), type = "character", help = "[Required] Output file name in .bed")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

df = read.csv(opts$csv)
df_out = df %>% filter(adj_beta_qval < 0.05) %>% select(c('phe_chr', 'phe_from', 'phe_to'))
write.table(df_out, file = opts$out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)