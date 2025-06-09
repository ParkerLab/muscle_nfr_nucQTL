#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)
library(dplyr)
library(stringr)

option_list <- list(
    make_option(c("--celltype"), type = "character", help = "[Required] The celltype that we are checking"),
    make_option(c("--overlap"), type = "character", help = "[Required] eQTL overlap info saved location"),
    make_option(c("--output"), type = "character", help = "[Required] eQTL overlap fisher exact test output filename")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

celltype = opts$celltype

df = read.table(opts$overlap)

# Create contingency table
contingency_table <- table(df$eqtl_overlap, df$nuc_coloc)
print(contingency_table)

# Run Fisher's Exact Test
result <- fisher.test(contingency_table)

# Print the result
print(result)

# save result
save(result, file = opts$output)