#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

option_list <- list(
  make_option(c("--input"), type = "character", help = "Row-wise Steiger TSV"),
  make_option(c("--celltype"), type = "character", help = "Cell type label for plot title"),
  make_option(c("--hits"), type = "character", default = "both", help = "Which hit set to plot: both, hit1, or hit2"),
  make_option(c("--out"), type = "character", help = "Output PNG"),
  make_option(c("--summary"), type = "character", default = NULL, help = "Optional output summary TSV")
)

opts <- parse_args(OptionParser(option_list = option_list))

if (is.null(opts$input) || is.null(opts$out) || is.null(opts$celltype)) {
  stop("--input, --celltype, and --out are required")
}

df <- read.table(opts$input, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

hit_cols <- switch(
  opts$hits,
  both = c("direction_hit1", "direction_hit2"),
  hit1 = c("direction_hit1"),
  hit2 = c("direction_hit2"),
  stop("--hits must be one of: both, hit1, hit2")
)

plot_df <- df %>%
  select(all_of(hit_cols)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "hit_set",
    values_to = "direction"
  ) %>%
  mutate(
    hit_set = recode(hit_set, direction_hit1 = "hit1", direction_hit2 = "hit2"),
    hit_set = factor(hit_set, levels = c("hit1", "hit2"))
  ) %>%
  filter(!is.na(direction), direction != "undetermined")

summary_df <- plot_df %>%
  count(hit_set, direction, name = "count") %>%
  arrange(hit_set, direction)

if (!is.null(opts$summary)) {
  write.table(summary_df, file = opts$summary, sep = "\t", quote = FALSE, row.names = FALSE)
}

if (nrow(summary_df) == 0) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No directionality calls available", size = 4) +
    theme_void() +
    labs(title = paste0(opts$celltype, " Steiger Directionality (", opts$hits, ")"))
  ggsave(filename = opts$out, plot = p, dpi = 300, width = 4, height = 4)
  quit(save = "no", status = 0)
}

direction_levels <- c("nfr_to_nuc", "nuc_to_nfr")
present_levels <- direction_levels[direction_levels %in% unique(summary_df$direction)]
if (length(setdiff(unique(summary_df$direction), direction_levels)) > 0) {
  present_levels <- c(present_levels, setdiff(unique(summary_df$direction), direction_levels))
}

summary_df$direction <- factor(summary_df$direction, levels = present_levels)

plot_colors <- c(
  nfr_to_nuc = "#FAA82D",
  nuc_to_nfr = "#455DAA"
)
plot_colors <- plot_colors[names(plot_colors) %in% present_levels]

p <- ggplot(summary_df, aes(x = hit_set, y = count, fill = direction)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(
    aes(label = count),
    position = position_dodge(width = 0.9),
    vjust = -0.15,
    size = 3
  ) +
  scale_fill_manual(values = plot_colors, drop = FALSE) +
  labs(
    x = "Lead SNP set",
    y = "Colocalized peak pairs",
    title = paste0(opts$celltype, " Steiger Directionality (", opts$hits, ")")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    plot.title = element_text(size = 10, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )

ggsave(filename = opts$out, plot = p, dpi = 300, width = 4, height = 4)
