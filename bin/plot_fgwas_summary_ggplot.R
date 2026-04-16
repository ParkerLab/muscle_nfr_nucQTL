#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(stringr)
  library(optparse)
})

option_list <- list(
  make_option(c("--input"), type = "character", dest = "input",
              help = "Path to wide fgwas summary TSV"),
  make_option(c("--outdir"), type = "character", dest = "outdir",
              help = "Output directory for figures"),
  make_option(c("--width"), type = "double", dest = "width", default = 10,
              help = "Figure width in inches [default %default]"),
  make_option(c("--height"), type = "double", dest = "height", default = 7,
              help = "Figure height in inches [default %default]")
)

opts <- parse_args(OptionParser(option_list = option_list))

if (is.null(opts$input) || opts$input == "") {
  stop("Missing required argument: --input", call. = FALSE)
}
if (is.null(opts$outdir) || opts$outdir == "") {
  stop("Missing required argument: --outdir", call. = FALSE)
}

dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

annotation_colors <- c(
  caPeaks = "#D55E00",
  nfrPeaks = "#E69F00",
  nucPeaks = "#0072B2"
)

check_required_columns <- function(df, columns, label) {
  missing_cols <- setdiff(columns, names(df))
  if (length(missing_cols) > 0) {
    stop(
      label, " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

save_plot_pair <- function(plot_obj, stem, outdir, width, height) {
  ggsave(
    filename = file.path(outdir, paste0(stem, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
  ggsave(
    filename = file.path(outdir, paste0(stem, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}

pivot_fgwas_coefficients <- function(df, annotations) {
  long_list <- lapply(annotations, function(annotation) {
    estimate_col <- paste0(annotation, "_ln")
    ci_low_col <- paste0(annotation, "_ln_lo")
    ci_high_col <- paste0(annotation, "_ln_hi")

    check_required_columns(
      df,
      c("trait", "model", estimate_col, ci_low_col, ci_high_col),
      paste0("Input table for ", annotation)
    )

    df %>%
      transmute(
        trait = .data$trait,
        model = .data$model,
        annotation = annotation,
        estimate = .data[[estimate_col]],
        ci_low = .data[[ci_low_col]],
        ci_high = .data[[ci_high_col]]
      )
  })

  bind_rows(long_list)
}

build_forest_plot <- function(df, annotations, title) {
  dodge <- position_dodge(width = 0.6)

  ggplot(df, aes(x = estimate, y = trait, color = annotation)) +
    geom_vline(xintercept = 0, color = "grey50", linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      position = dodge,
      width = 0.2,
      linewidth = 0.6,
      orientation = "y"
    ) +
    geom_point(
      position = dodge,
      size = 2.4
    ) +
    scale_color_manual(values = annotation_colors[annotations], drop = FALSE) +
    labs(
      x = "fgwas log-enrichment coefficient",
      y = "Trait",
      color = "Annotation",
      title = title
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

fgwas_df <- read_tsv(opts$input, show_col_types = FALSE)

check_required_columns(
  fgwas_df,
  c(
    "trait", "model",
    "caPeaks_ln", "caPeaks_ln_lo", "caPeaks_ln_hi",
    "nfrPeaks_ln", "nfrPeaks_ln_lo", "nfrPeaks_ln_hi",
    "nucPeaks_ln", "nucPeaks_ln_lo", "nucPeaks_ln_hi"
  ),
  "Input table"
)

trait_order <- fgwas_df %>%
  filter(model == "ca_only") %>%
  arrange(desc(caPeaks_ln)) %>%
  pull(trait)

if (length(trait_order) == 0) {
  trait_order <- fgwas_df %>%
    distinct(trait) %>%
    pull(trait)
}

coeff_long <- pivot_fgwas_coefficients(fgwas_df, c("caPeaks", "nfrPeaks", "nucPeaks")) %>%
  mutate(
    trait = factor(trait, levels = rev(unique(trait_order))),
    annotation = factor(annotation, levels = c("caPeaks", "nfrPeaks", "nucPeaks"))
  )

marginal_df <- bind_rows(
  coeff_long %>% filter(model == "ca_only", annotation == "caPeaks"),
  coeff_long %>% filter(model == "nfr_only", annotation == "nfrPeaks"),
  coeff_long %>% filter(model == "nuc_only", annotation == "nucPeaks")
) %>%
  filter(!is.na(estimate), !is.na(ci_low), !is.na(ci_high))

ca_nfr_df <- coeff_long %>%
  filter(model == "ca_nfr", annotation %in% c("caPeaks", "nfrPeaks")) %>%
  filter(!is.na(estimate), !is.na(ci_low), !is.na(ci_high)) %>%
  mutate(annotation = fct_drop(annotation))

marginal_plot <- build_forest_plot(
  marginal_df,
  annotations = c("caPeaks", "nfrPeaks", "nucPeaks"),
  title = "Marginal fgwas Log-Enrichment Coefficients"
)

ca_nfr_plot <- build_forest_plot(
  ca_nfr_df,
  annotations = c("caPeaks", "nfrPeaks"),
  title = "Conditional fgwas Log-Enrichment Coefficients: ca_nfr"
)

print(marginal_plot)
print(ca_nfr_plot)

save_plot_pair(marginal_plot, "fgwas_forest_marginal", opts$outdir, opts$width, opts$height)
save_plot_pair(ca_nfr_plot, "fgwas_forest_ca_nfr", opts$outdir, opts$width, opts$height)
