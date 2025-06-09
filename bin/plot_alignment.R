#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(coloc)
library(dplyr)
library(optparse)

option_list <- list(
    make_option(c("--filepath"), type = "character", help = "[Required] Filepath for the coloc dataset"),
    make_option(c("--out_filepath"), type = "character", help = "[Required] Which directory to save the output plots")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# variables
filepath = opts$filepath
out_filepath = opts$out_filepath

load(file=filepath)

bprod=outer(coloc_d1$beta/sqrt(coloc_d1$varbeta),coloc_d1$beta/sqrt(coloc_d1$varbeta),"*")
## plot(bprod,D$LD[D$snp,D$snp],xlab="product of z scores",ylab="LD")
tmp=(bprod/coloc_d1$LD)[abs(coloc_d1$LD) > 0.2]
png(paste0(out_filepath, ".png"))
hist(tmp,
   xlab="ratio of product of Z scores to LD \n",
   main="alignment check plot",
   sub="\n expect most values to be positive \n symmetry is a warning sign of potentially poor alignment")
legend("topright",legend=paste0("% positive = ",100*round(mean( tmp > 0 ),3)))
abline(v=0,col="red")
dev.off()