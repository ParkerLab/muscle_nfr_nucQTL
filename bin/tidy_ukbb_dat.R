#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(susieR)
library(dplyr)
library(vcfR)
library(coloc)
library(Rfast)
library(optparse)

option_list <- list(
    make_option(c("--qtl_file"), type = "character", help = "[Required] Lead variant output filename"),
    make_option(c("--ukbb_path"), type = "character", help = "[Required] Using which VCF information, ukbb?"),
    make_option(c("--min_corr"), type = "numeric", help = "[Required] min abs corr between any snp pairs"),
    make_option(c("--num_L"), type = "numeric", help = "[Required] Max number of LD vars"),
    make_option(c("--max_it"), type = "numeric", help = "[Required] Max number of susie iterations"),
    make_option(c("--outdir"), type = "character", help = "[Required] Susie models output directory")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

