#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(dplyr)
library(stringr)
library(optparse)

option_list <- list(
    make_option(c("--method"), type = "character", help = "[Required] Correcting phenotypes from which method"),
    make_option(c("--outdir"), type = "character", help = "[Required] output directory")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# input variables
# the method
method = opts$method
# outputs dir
outdir = opts$outdir

# -----------process inputs-----------
annotation = strsplit(method, split = "_counts")[[1]][1]

bed_file = paste0(strsplit(method, split = "_F")[[1]][1], ".bed.gz")
print(bed_file)
command_correct = "QTLtools correct --bed ../../data/QTL_phenotype_data/"
command_correct = paste0(command_correct, bed_file)
command_correct = paste0(command_correct, " --out ")
command_correct = paste0(command_correct, outdir)
command_correct = paste0(command_correct, annotation)
command_correct = paste0(command_correct, "_corrected.bed.gz --cov ../../data/QTL_covariates_data/")
command_correct = paste0(command_correct, method)
command_correct = paste0(command_correct, "_cov.txt")
print(command_correct)
system(command_correct, intern=TRUE)

# read in and process the corrected bed file
bed_corrected = paste0(outdir, annotation)
bed_corrected_bed = paste0(bed_corrected, "_corrected.bed.gz")
dat = as.data.frame(read.table(bed_corrected_bed, comment.char=''), header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
# assgin headers
names(dat) <- dat[1,]
dat <- dat[-1,]
# save
print(paste0(bed_corrected, "_corrected.Rda"))
save(dat, file = paste0(bed_corrected, "_corrected.Rda"))