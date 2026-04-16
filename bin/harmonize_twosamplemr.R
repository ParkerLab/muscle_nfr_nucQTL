#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(TwoSampleMR)
})

option_list <- list(
  make_option(c("--nfr"), type = "character", help = "Annotated nfr exposure table"),
  make_option(c("--nuc"), type = "character", help = "Annotated nuc outcome table"),
  make_option(c("--n_exp"), type = "integer", help = "Exposure sample size"),
  make_option(c("--n_out"), type = "integer", help = "Outcome sample size"),
  make_option(c("--out"), type = "character", help = "Output harmonized TSV")
)

opts <- parse_args(OptionParser(option_list = option_list))

exp <- read.table(opts$nfr, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
out <- read.table(opts$nuc, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

req_exp <- c("hit1", "slope", "slope_se", "nom_pval", "effect_allele", "other_allele", "eaf")
req_out <- c("hit2", "slope", "slope_se", "nom_pval", "effect_allele", "other_allele", "eaf")

miss_exp <- setdiff(req_exp, colnames(exp))
miss_out <- setdiff(req_out, colnames(out))
if (length(miss_exp) > 0) stop(paste("Missing exposure cols:", paste(miss_exp, collapse = ",")))
if (length(miss_out) > 0) stop(paste("Missing outcome cols:", paste(miss_out, collapse = ",")))

exp2 <- data.frame(
  SNP = exp$hit1,
  beta = as.numeric(exp$slope),
  se = as.numeric(exp$slope_se),
  pval = as.numeric(exp$nom_pval),
  effect_allele = exp$effect_allele,
  other_allele = exp$other_allele,
  eaf = as.numeric(exp$eaf),
  samplesize = as.integer(opts$n_exp),
  id.exposure = "nfr_splitA",
  exposure = "nfr_splitA"
)

out2 <- data.frame(
  SNP = out$hit2,
  beta = as.numeric(out$slope),
  se = as.numeric(out$slope_se),
  pval = as.numeric(out$nom_pval),
  effect_allele = out$effect_allele,
  other_allele = out$other_allele,
  eaf = as.numeric(out$eaf),
  samplesize = as.integer(opts$n_out),
  id.outcome = "nuc_splitB",
  outcome = "nuc_splitB"
)

exp2 <- exp2[!duplicated(exp2$SNP), ]
out2 <- out2[!duplicated(out2$SNP), ]

exp_dat <- format_data(
  exp2,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  samplesize_col = "samplesize"
)

out_dat <- format_data(
  out2,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  samplesize_col = "samplesize"
)

harm <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat,
  action = 2
)

write.table(harm, file = opts$out, sep = "\t", quote = FALSE, row.names = FALSE)
