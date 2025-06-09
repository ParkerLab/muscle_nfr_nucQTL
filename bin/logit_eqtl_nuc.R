#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(dplyr)
library(ggplot2)
library(optparse)

option_list <- list(
    make_option(c("--celltype"), type = "character", help = "[Required] Looking at which celltype"),
    make_option(c("--input"), type = "character", help = "[Required] input file with eqtl overlap results"),
    make_option(c("--nfr_susie_dir"), type = "character", help = "[Required] directory of nfr susie results"),
    make_option(c("--outdir"), type = "character", help = "[Required] output directory")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

celltype = opts$celltype

dir = opts$input
# read in the eQTL overlap results
df = read.table(dir)
outdir_nfr = paste0(opts$nfr_susie_dir, celltype, "_nfr_susie/")
# Initialize a new column in df to store the mean values
df$mean_nfr <- NA  # Using NA as a placeholder for missing values

# Loop over each row of df to calculate the mean for the specified NFR_peak
for (i in 1:nrow(df)) {
  # Extract the NFR_peak identifier from the current row of df
    pheid <- df[i, ]$NFR_peak
    chr =df[i, ]$Chr
    load(file = paste0(outdir_nfr, chr, '/', pheid, '_susie_dat.Rda'))
    # add NFR y
    df[i,]$mean_nfr = mean(coloc_d1$y)
}

# Run logistic regression
logistic_model <- glm(eqtl_overlap ~ nuc_coloc + mean_nfr, 
                      data = df, 
                      family = binomial(link = "logit"))

# Display the summary of the model to view results
summary(logistic_model)
# Save the entire model object to a file
saveRDS(logistic_model, file = opts$outdir)