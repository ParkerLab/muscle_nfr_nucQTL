#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)
library(dplyr)

option_list <- list(
    make_option(c("--celltype"), type = "character", help = "[Required] Which celltype are we looking at"),
    make_option(c("--region"), type = "character", help = "[Required] Which region? nfr or nuc?"),
    make_option(c("--input_count"), type = "character", help = "[Required] Where is this raw count matrix"),
    make_option(c("--output"), type = "character", help = "[Required] Where to store the tidied count matrix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# input
celltype = opts$celltype
region = opts$region
input_count = opts$input_count
output = opts$output

# read in featureCounts output file
d_ = read.table(input_count, header=TRUE)
# remove sex chromosome
d_ = d_ %>% filter(Chr != 'chrY') %>% filter(Chr != 'chrX')
# extract just counts columns
sample_cols = names(d_)[ grep(pattern=".bam$", x=names(d_), perl=TRUE)]
d = d_[,sample_cols]
# add peak identifiers
row.names(d) = d_$Geneid
# reformat sample identifiers
## identify pattern
pattern = paste0('X.home.xiaoouw.Muscle_snATAC_nucQTL.data.bams.by.cluster.sample.', 
                 celltype, 
                 #'.',
                 #region,
                 '.atac.',
                 celltype,
                 '.fusion.')
## rename
names(d) <- gsub(".", "_", sapply(strsplit(sample_cols, pattern), "[[", 2), fixed=TRUE)
# convert to matrix
m = as.matrix(d)
#save as rds
saveRDS(m, file=output)