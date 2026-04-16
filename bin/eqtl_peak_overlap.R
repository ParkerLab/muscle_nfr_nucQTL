#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)

option_list <- list(
    make_option(c("--celltype"), type = "character", help = "[Required] The celltype that we are checking"),
    make_option(c("--eqtl"), type = "character", help = "[Required] eQTL summary table from eqtl_summary.R"),
    make_option(c("--qtl_bed"), type = "character", help = "[Required] target QTL peak BED file used to resolve SuSiE credible sets"),
    make_option(c("--qtl_susie_root"), type = "character", help = "[Required] root directory containing per-peak SuSiE credible sets"),
    make_option(c("--qtl_label"), type = "character", help = "[Required] label for the target QTL class"),
    make_option(c("--output"), type = "character", help = "[Required] output summary filename")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

read_peak_bed <- function(path){
    if (!file.exists(path)){
        stop(paste0("Missing input file: ", path))
    }
    read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
}

peak_name_candidates <- function(chr, start, end){
    start_num <- suppressWarnings(as.integer(start))
    end_num <- suppressWarnings(as.integer(end))
    candidates <- c()

    if (!is.na(start_num) && !is.na(end_num)){
        candidates <- c(
            paste0(chr, ".", start_num - 1L, ".", end_num),
            paste0(chr, ".", start_num, ".", end_num)
        )
    } else {
        candidates <- paste0(chr, ".", start, ".", end)
    }

    unique(candidates)
}

extract_eqtl_ids <- function(eqtl_bed){
    eqtl_df <- read.table(eqtl_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
    if (nrow(eqtl_df) == 0){
        return(character())
    }

    ids <- unique(c(
        as.character(eqtl_df[[4]]),
        as.character(eqtl_df[[7]])
    ))
    ids <- ids[!is.na(ids) & ids != ""]
    unique(ids)
}

resolve_qtl_ids <- function(qtl_bed, qtl_susie_root){
    peak_df <- read_peak_bed(qtl_bed)
    if (nrow(peak_df) == 0){
        stop(paste0("No peaks found in ", qtl_bed))
    }

    qtl_ids <- character()
    for (i in seq_len(nrow(peak_df))){
        chr <- as.character(peak_df[i, 1])
        start <- as.character(peak_df[i, 2])
        end <- as.character(peak_df[i, 3])
        cset_path <- NULL
        peak_cset <- NULL

        for (peak_name in peak_name_candidates(chr, start, end)){
            candidate_path <- file.path(qtl_susie_root, chr, paste0(peak_name, ".susie_cset95.tsv"))
            if (file.exists(candidate_path)){
                cset_path <- candidate_path
                peak_cset <- read.table(candidate_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
                break
            }
        }

        if (is.null(cset_path)){
            warning(paste0("Skipping missing qtl credible set: ", file.path(qtl_susie_root, chr, paste0(peak_name_candidates(chr, start, end)[1], ".susie_cset95.tsv"))))
            next
        }

        if (nrow(peak_cset) == 0 || !"snp" %in% colnames(peak_cset)){
            next
        }

        qtl_ids <- c(qtl_ids, as.character(peak_cset$snp))
    }

    qtl_ids <- qtl_ids[!is.na(qtl_ids) & qtl_ids != ""]
    unique(qtl_ids)
}

resolve_qtl_ids_from_root <- function(qtl_susie_root){
    chr_dirs <- list.dirs(qtl_susie_root, recursive = FALSE, full.names = TRUE)
    chr_dirs <- chr_dirs[grepl("/chr[0-9XY]+$", chr_dirs)]

    qtl_ids <- character()
    for (chr_dir in chr_dirs){
        chr <- basename(chr_dir)
        summary_path <- file.path(chr_dir, paste0(chr, "_susie_cset95_summary.tsv"))
        if (!file.exists(summary_path)){
            next
        }

        summary_df <- read.table(summary_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
        if (nrow(summary_df) == 0){
            next
        }

        for (i in seq_len(nrow(summary_df))){
            peak_chr <- as.character(summary_df[i, ]$Chr)
            peak_start <- as.character(summary_df[i, ]$Start)
            peak_end <- as.character(summary_df[i, ]$End)

            if (summary_df[i, ]$n95cset == 0){
                next
            }

            peak_name_list <- peak_name_candidates(peak_chr, peak_start, peak_end)
            if (summary_df[i, ]$n95cset == 1){
                cset_path <- NULL
                peak_cset <- NULL
                for (peak_name in peak_name_list){
                    candidate_path <- file.path(chr_dir, paste0(peak_name, ".susie_cset95.tsv"))
                    if (file.exists(candidate_path)){
                        cset_path <- candidate_path
                        peak_cset <- read.table(candidate_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
                        break
                    }
                }
                if (is.null(cset_path)){
                    warning(paste0("Skipping missing qtl credible set: ", file.path(chr_dir, paste0(peak_name_list[1], ".susie_cset95.tsv"))))
                    next
                }
                if (nrow(peak_cset) > 0 && "snp" %in% colnames(peak_cset)){
                    qtl_ids <- c(qtl_ids, as.character(peak_cset$snp))
                }
            } else {
                susie_path <- NULL
                pip_path <- NULL
                peak_name <- NULL
                for (candidate_name in peak_name_list){
                    candidate_susie <- file.path(chr_dir, paste0(candidate_name, ".susie.rda"))
                    candidate_pip <- file.path(chr_dir, paste0(candidate_name, ".susie_pipdf.tsv"))
                    if (file.exists(candidate_susie) && file.exists(candidate_pip)){
                        peak_name <- candidate_name
                        susie_path <- candidate_susie
                        pip_path <- candidate_pip
                        break
                    }
                }
                if (is.null(peak_name)){
                    warning(paste0("Skipping multi-set peak with missing files: ", peak_name_list[1]))
                    next
                }

                load(susie_path)
                pip_df <- read.table(pip_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
                snp_list <- c()
                for (j in seq_along(S1$sets$cs)){
                    snp_list <- c(snp_list, S1$sets$cs[[j]])
                }
                snp_df <- pip_df[snp_list, , drop = FALSE]
                if (nrow(snp_df) > 0 && "snp" %in% colnames(snp_df)){
                    qtl_ids <- c(qtl_ids, as.character(snp_df$snp))
                }
            }
        }
    }

    qtl_ids <- qtl_ids[!is.na(qtl_ids) & qtl_ids != ""]
    unique(qtl_ids)
}

eqtl_summary_df <- read.table(opts$eqtl, sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")
required_cols <- c("gene", "chr", "pos", "leadvar", "path")
missing_cols <- setdiff(required_cols, colnames(eqtl_summary_df))
if (length(missing_cols) > 0){
    stop(paste0("Missing columns in eQTL summary table: ", paste(missing_cols, collapse = ", ")))
}

qtl_ids <- unique(c(
    resolve_qtl_ids(opts$qtl_bed, opts$qtl_susie_root),
    resolve_qtl_ids_from_root(opts$qtl_susie_root)
))
if (length(qtl_ids) == 0){
    stop(paste0("No qtl identifiers could be resolved from ", opts$qtl_bed))
}

overlap_vec <- logical(nrow(eqtl_summary_df))
for (i in seq_len(nrow(eqtl_summary_df))){
    eqtl_ids <- extract_eqtl_ids(eqtl_summary_df$path[i])
    overlap_vec[i] <- length(intersect(eqtl_ids, qtl_ids)) > 0
}

overlap_count <- sum(overlap_vec)
total_eqtls <- nrow(eqtl_summary_df)
overlap_percent <- if (total_eqtls == 0) NA_real_ else (overlap_count / total_eqtls) * 100

summary_df <- data.frame(
    celltype = opts$celltype,
    qtl_label = opts$qtl_label,
    total_eqtls = total_eqtls,
    overlapping_eqtls = overlap_count,
    overlap_percent = overlap_percent,
    qtl_bed = opts$qtl_bed,
    stringsAsFactors = FALSE
)

write.table(summary_df, file = opts$output, sep = "\t", quote = FALSE, row.names = FALSE)
