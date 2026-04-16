#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

library(optparse)
library(dplyr)
library(ggplot2)
library(patchwork)

option_list <- list(
    make_option(c("--gwas_peak"), type = "character", help = "[Required] GWAS peak identifier, e.g. chr1-100-101"),
    make_option(c("--celltype"), type = "character", default = "Type_2a", help = "[Optional] Cell type label"),
    make_option(c("--gwas_trait"), type = "character", default = "grip_strength", help = "[Optional] GWAS trait identifier"),
    make_option(c("--gwas_label"), type = "character", default = "grip_strength", help = "[Optional] GWAS label"),
    make_option(c("--coloc_nfr"), type = "character", help = "[Required] nfr coloc summary Rda"),
    make_option(c("--coloc_nuc"), type = "character", help = "[Required] nuc coloc summary Rda"),
    make_option(c("--coloc_ca_np"), type = "character", help = "[Required] ca_np coloc summary Rda"),
    make_option(c("--nfr_susie_root"), type = "character", help = "[Required] Root directory for nfr susie_dat files"),
    make_option(c("--nuc_susie_root"), type = "character", help = "[Required] Root directory for nuc susie_dat files"),
    make_option(c("--ca_susie_root"), type = "character", help = "[Required] Root directory for ca_np-as-ca susie_dat files"),
    make_option(c("--gwas_susie_root"), type = "character", help = "[Required] Root directory for GWAS susie_dat files"),
    make_option(c("--out_png"), type = "character", help = "[Required] Output PNG path")
)

opts <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

base_text_size <- 0.7
base_point_size <- 0.55

read_rda_object <- function(filepath, object_pattern = NULL) {
    env <- new.env(parent = emptyenv())
    load(filepath, envir = env)
    object_names <- ls(env)

    if (!is.null(object_pattern)) {
        matched <- grep(object_pattern, object_names, value = TRUE)
        if (length(matched) > 0) {
            return(get(matched[[1]], envir = env))
        }
    }

    if (length(object_names) != 1) {
        stop("Could not unambiguously resolve object in ", filepath)
    }

    get(object_names[[1]], envir = env)
}

normalize_summary_df <- function(df) {
    if (is.null(df) || nrow(df) == 0) {
        return(data.frame())
    }

    if ("PP.H4.abf" %in% names(df)) {
        df <- df %>% arrange(desc(.data$`PP.H4.abf`))
    }

    df
}

pick_best_signal <- function(filepath) {
    df <- read_rda_object(filepath, "^df_res")
    df <- normalize_summary_df(df)
    if (nrow(df) == 0) {
        return(NULL)
    }
    df[1, , drop = FALSE]
}

find_peak_column <- function(df) {
    peak_candidates <- c("QTL_peak", "NFR_peak", "Nucleosomal_peak", "peak_id")
    matched <- peak_candidates[peak_candidates %in% names(df)]
    if (length(matched) == 0) {
        stop("Could not identify QTL peak column from: ", paste(names(df), collapse = ", "))
    }
    matched[[1]]
}

peak_chr <- function(peak_id) {
    if (grepl("\\.", peak_id)) {
        return(strsplit(peak_id, "\\.")[[1]][1])
    }
    strsplit(peak_id, "-")[[1]][1]
}

load_susie_dat <- function(filepath) {
    dat <- read_rda_object(filepath, "^coloc_d1$")
    required <- c("snp", "pvalues", "position", "LD")
    missing <- setdiff(required, names(dat))
    if (length(missing) > 0) {
        stop("Missing fields in ", filepath, ": ", paste(missing, collapse = ", "))
    }
    dat
}

resolve_ld_column <- function(ld_matrix, snp_id) {
    ld_df <- as.data.frame(ld_matrix)
    if (!snp_id %in% colnames(ld_df)) {
        return(NULL)
    }
    ld_df[[snp_id]]
}

build_track_df <- function(filepath, ld_reference_snp, fallback_snp, track_label) {
    dat <- load_susie_dat(filepath)
    ld_values <- resolve_ld_column(dat$LD, ld_reference_snp)
    if (is.null(ld_values)) {
        ld_values <- rep(NA_real_, length(dat$snp))
    }

    tibble(
        rsid = dat$snp,
        pval = dat$pvalues,
        position = dat$position,
        r2 = suppressWarnings(as.numeric(ld_values)^2),
        lead = dat$snp == fallback_snp,
        ld_reference = ld_reference_snp,
        track = track_label
    )
}

empty_track_plot <- function(track_label, subtitle = "No coloc peak found") {
    ggplot() +
        annotate("text", x = 0, y = 0, label = subtitle, size = 4 * base_text_size) +
        labs(title = track_label, x = NULL, y = "-log10(P)") +
        theme_void() +
        theme(
            plot.title = element_text(hjust = 0, face = "bold", size = 11 * base_text_size),
            plot.margin = margin(5.5, 10, 5.5, 10)
    )
}

track_plot <- function(df, track_label, y_label, x_limits, highlight_snp, subtitle = NULL, show_x = FALSE) {
    df <- df %>%
        mutate(
            neg_log10_p = -log10(pval),
            highlight = rsid == highlight_snp
        )

    lead_df <- df %>% filter(highlight)
    p <- ggplot(df, aes(x = position, y = neg_log10_p)) +
        geom_point(aes(color = r2), size = 1.8 * base_point_size, na.rm = TRUE) +
        geom_point(
            data = df %>% filter(is.na(r2)),
            color = "grey70",
            size = 1.6 * base_point_size,
            inherit.aes = FALSE,
            aes(x = position, y = neg_log10_p)
        ) +
        scale_color_gradient(
            low = "#d9d9d9",
            high = "#08519c",
            limits = c(0, 1),
            na.value = "grey70"
        ) +
        coord_cartesian(xlim = x_limits) +
        labs(
            title = track_label,
            subtitle = subtitle,
            x = if (show_x) "Genomic position (bp)" else NULL,
            y = y_label,
            color = expression(r^2)
        ) +
        guides(
            color = guide_colorbar(
                barwidth = unit(3, "mm"),
                barheight = unit(18, "mm"),
                ticks = FALSE
            )
        ) +
        theme_classic() +
        theme(
            plot.title = element_text(face = "bold", size = 11 * base_text_size),
            plot.subtitle = element_text(size = 9 * base_text_size),
            axis.title.x = element_text(size = 10 * base_text_size),
            axis.title.y = element_text(size = 10 * base_text_size),
            axis.text.x = element_text(size = 7 * base_text_size),
            axis.text.y = element_text(size = 7 * base_text_size),
            legend.title = element_text(size = 9 * base_text_size),
            legend.text = element_text(size = 7 * base_text_size),
            legend.key.height = unit(18, "mm"),
            legend.key.width = unit(3, "mm"),
            legend.position = "right",
            plot.margin = margin(5.5, 10, 5.5, 10)
        )

    if (nrow(lead_df) > 0) {
        p <- p +
            geom_point(data = lead_df, color = "#cb181d", size = 2.2 * base_point_size) +
            geom_text(
                data = lead_df,
                aes(label = rsid),
                nudge_y = 0.35,
                size = 2.8 * base_text_size,
                color = "#cb181d",
                check_overlap = TRUE
            )
    }

    if (!show_x) {
        p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }

    p
}

resolve_qtl_filepath <- function(root_dir, peak_id) {
    file.path(root_dir, peak_chr(peak_id), paste0(peak_id, "_susie_dat.Rda"))
}

nfr_signal <- pick_best_signal(opts$coloc_nfr)
nuc_signal <- pick_best_signal(opts$coloc_nuc)
ca_signal <- pick_best_signal(opts$coloc_ca_np)

reference_signal <- dplyr::bind_rows(nfr_signal, nuc_signal, ca_signal) %>%
    arrange(desc(.data$`PP.H4.abf`)) %>%
    slice_head(n = 1)

if (nrow(reference_signal) == 0) {
    stop("No coloc signals found for ", opts$gwas_peak)
}

gwas_lead_snp <- reference_signal$hit1[[1]]

build_qtl_track <- function(signal_df, root_dir, track_label, y_label) {
    if (is.null(signal_df) || nrow(signal_df) == 0) {
        return(list(plot = empty_track_plot(track_label), data = NULL))
    }

    peak_col <- find_peak_column(signal_df)
    peak_id <- signal_df[[peak_col]][[1]]
    qtl_lead_snp <- signal_df$hit2[[1]]
    pph4 <- round(signal_df$`PP.H4.abf`[[1]], 4)
    subtitle <- paste0(peak_id, " | PPH4=", format(pph4, nsmall = 4))
    filepath <- resolve_qtl_filepath(root_dir, peak_id)

    if (!file.exists(filepath)) {
        return(list(
            plot = empty_track_plot(track_label, paste("Missing", basename(filepath))),
            data = NULL
        ))
    }

    df <- build_track_df(filepath, gwas_lead_snp, qtl_lead_snp, track_label)
    plot <- track_plot(
        df = df,
        track_label = track_label,
        y_label = y_label,
        x_limits = c(0, 1),
        highlight_snp = qtl_lead_snp,
        subtitle = paste0(subtitle, " | LD ref: ", gwas_lead_snp),
        show_x = FALSE
    )

    list(plot = plot, data = df)
}

nfr_track <- build_qtl_track(nfr_signal, opts$nfr_susie_root, "nfrQTL", "nfrQTL\n-log10(P)")
nuc_track <- build_qtl_track(nuc_signal, opts$nuc_susie_root, "nucQTL", "nucQTL\n-log10(P)")
ca_track <- build_qtl_track(ca_signal, opts$ca_susie_root, "caQTL", "caQTL\n-log10(P)")

gwas_filepath <- file.path(opts$gwas_susie_root, paste0("gwas_", opts$gwas_peak, "_susie_dat.Rda"))
if (!file.exists(gwas_filepath)) {
    stop("Missing GWAS susie_dat file: ", gwas_filepath)
}
gwas_track_df <- build_track_df(gwas_filepath, gwas_lead_snp, gwas_lead_snp, opts$gwas_label)

track_data <- list(nfr_track$data, nuc_track$data, ca_track$data, gwas_track_df)
track_data <- track_data[!vapply(track_data, is.null, logical(1))]
x_positions <- unlist(lapply(track_data, function(x) x$position))
x_limits <- range(x_positions, na.rm = TRUE)

nfr_plot <- if (!is.null(nfr_track$data)) {
    nfr_track$plot + coord_cartesian(xlim = x_limits)
} else {
    nfr_track$plot
}
nuc_plot <- if (!is.null(nuc_track$data)) {
    nuc_track$plot + coord_cartesian(xlim = x_limits)
} else {
    nuc_track$plot
}
ca_plot <- if (!is.null(ca_track$data)) {
    ca_track$plot + coord_cartesian(xlim = x_limits)
} else {
    ca_track$plot
}
gwas_plot <- track_plot(
    df = gwas_track_df,
    track_label = opts$gwas_label,
    y_label = "GWAS\n-log10(P)",
    x_limits = x_limits,
    highlight_snp = gwas_lead_snp,
    subtitle = paste0("lead SNP: ", gwas_lead_snp, " | LD ref: ", gwas_lead_snp),
    show_x = TRUE
)

combined_plot <- nfr_plot / nuc_plot / ca_plot / gwas_plot +
    plot_annotation(
        title = paste(opts$celltype, opts$gwas_trait, opts$gwas_peak, sep = " | "),
        theme = theme(
            plot.title = element_text(size = 12 * base_text_size, face = "bold")
        )
    )

ggsave(
    filename = opts$out_png,
    plot = combined_plot,
    width = 4,
    height = 5,
    dpi = 300
)
