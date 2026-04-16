#!/usr/bin/env Rscript

options(stringsAsFactors=FALSE)
library(optparse)

option_list <- list(
    make_option(c("--inputs"), type = "character", help = "[Required] Space-separated list of summary TSV files"),
    make_option(c("--output"), type = "character", help = "[Required] output combined summary filename")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

input_files <- unlist(strsplit(opts$inputs, "\\s+"))
input_files <- input_files[input_files != ""]
input_files <- unique(input_files)

if (length(input_files) == 0){
    stop("No input summary files were provided")
}

summary_list <- list()
for (input_file in input_files){
    if (!file.exists(input_file)){
        stop(paste0("Missing input summary file: ", input_file))
    }

    summary_df <- read.table(input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    summary_df$summary_file <- basename(input_file)
    summary_list[[length(summary_list) + 1]] <- summary_df
}

combined_df <- do.call(rbind, summary_list)
write.table(combined_df, file = opts$output, sep = "\t", quote = FALSE, row.names = FALSE)
