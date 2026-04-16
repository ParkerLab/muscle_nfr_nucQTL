#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(TwoSampleMR)
})

option_list <- list(
  make_option(c("--harmonized"), type = "character", help = "Input harmonized TwoSampleMR table"),
  make_option(c("--pval_threshold"), type = "double", default = 5e-6, help = "Exposure p-value threshold"),
  make_option(c("--out"), type = "character", help = "Output Steiger results TSV")
)

opts <- parse_args(OptionParser(option_list = option_list))

df <- read.table(opts$harmonized, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

if (!"mr_keep" %in% colnames(df)) {
  stop("harmonized input missing mr_keep column")
}
if (!"pval.exposure" %in% colnames(df)) {
  stop("harmonized input missing pval.exposure column")
}

# First-pass instrument filtering: keep harmonized SNPs with strong exposure association.
df_f <- df[df$mr_keep & !is.na(df$pval.exposure) & df$pval.exposure <= opts$pval_threshold, ]

if (nrow(df_f) == 0) {
  out <- data.frame(
    n_input = nrow(df),
    n_after_filter = 0,
    pval_threshold = opts$pval_threshold,
    clumping_applied = FALSE,
    note = "No SNPs passed filter"
  )
  write.table(out, file = opts$out, sep = "\t", quote = FALSE, row.names = FALSE)
  quit(save = "no", status = 0)
}

# Steiger directionality test using harmonized exposure/outcome summary stats.
steiger <- steiger_filtering(df_f)

# Add run-level metadata for traceability.
steiger$input_n <- nrow(df)
steiger$filtered_n <- nrow(df_f)
steiger$pval_threshold <- opts$pval_threshold
steiger$clumping_applied <- FALSE

write.table(steiger, file = opts$out, sep = "\t", quote = FALSE, row.names = FALSE)
