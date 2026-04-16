#!/usr/bin/env Rscript

library(optparse)

option_list <- list(
  make_option(c("--norm_rds"), type = "character", help = "Input normalized matrix RDS"),
  make_option(c("--keep_peaks_file"), type = "character", help = "Peak list file (TSV/text)"),
  make_option(c("--keep_peaks_col"), type = "character", default = NULL, help = "Column in keep_peaks_file"),
  make_option(c("--out_rds"), type = "character", help = "Output subset RDS")
)

opts <- parse_args(OptionParser(option_list = option_list))

mat <- readRDS(opts$norm_rds)
if (is.null(dim(mat))) {
  stop("Normalized RDS is not a matrix/data.frame")
}

if (is.null(opts$keep_peaks_col)) {
  peaks <- readLines(opts$keep_peaks_file)
  peaks <- peaks[peaks != ""]
} else {
  kp <- read.table(opts$keep_peaks_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!(opts$keep_peaks_col %in% colnames(kp))) {
    stop(paste("keep_peaks_col not found:", opts$keep_peaks_col))
  }
  peaks <- as.character(kp[[opts$keep_peaks_col]])
  peaks <- peaks[!is.na(peaks) & peaks != ""]
}

peak_set <- unique(peaks)

# normalization.R output is sample x peak in this project; subset by columns.
available <- colnames(mat)
if (is.null(available)) {
  stop("Input matrix has no colnames; cannot subset by peak IDs")
}

keep <- intersect(peak_set, available)
missing <- setdiff(peak_set, available)

if (length(keep) == 0) {
  stop("No requested peaks found in normalized matrix")
}

mat_sub <- mat[, keep, drop = FALSE]

if (length(missing) > 0) {
  message(sprintf("WARNING: %d requested peaks missing from normalized matrix", length(missing)))
  message(paste("Missing (first 10):", paste(head(missing, 10), collapse = ", ")))
}

saveRDS(mat_sub, opts$out_rds)
message(sprintf("Wrote subset matrix with %d samples x %d peaks", nrow(mat_sub), ncol(mat_sub)))
