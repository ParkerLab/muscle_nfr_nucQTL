#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(TwoSampleMR)
  library(qvalue)
})

option_list <- list(
  make_option(c("--input"), type = "character", help = "Row-wise Steiger input TSV"),
  make_option(c("--q-threshold"), type = "double", default = 0.05, help = "Q-value threshold for direction assignment"),
  make_option(c("--out"), type = "character", help = "Output row-wise directionality TSV")
)

opts <- parse_args(OptionParser(option_list = option_list))
if (is.null(opts$input) || is.null(opts$out)) {
  stop("Both --input and --out are required")
}

df <- read.table(opts$input, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

num_cols <- c(
  "h1_eaf", "h2_eaf",
  "h1_nfr_beta", "h1_nfr_se", "h1_nfr_pval", "h1_nuc_beta", "h1_nuc_se", "h1_nuc_pval",
  "h2_nfr_beta", "h2_nfr_se", "h2_nfr_pval", "h2_nuc_beta", "h2_nuc_se", "h2_nuc_pval",
  "samplesize_nfr", "samplesize_nuc"
)
for (cn in intersect(num_cols, colnames(df))) {
  df[[cn]] <- suppressWarnings(as.numeric(df[[cn]]))
}

pick_first_col <- function(dat, candidates) {
  for (nm in candidates) {
    if (nm %in% colnames(dat)) return(dat[[nm]][1])
  }
  NA
}

run_directionality_one <- function(
  snp,
  ea,
  oa,
  eaf,
  beta_nfr,
  se_nfr,
  p_nfr,
  n_nfr,
  beta_nuc,
  se_nuc,
  p_nuc,
  n_nuc
) {
  needed <- c(beta_nfr, se_nfr, p_nfr, n_nfr, beta_nuc, se_nuc, p_nuc, n_nuc)
  if (is.na(ea) || is.na(oa) || any(is.na(needed))) {
    return(list(steiger_dir = NA, steiger_p = NA_real_, r2_nfr = NA_real_, r2_nuc = NA_real_, note = "missing"))
  }

  exp_df <- data.frame(
    SNP = as.character(snp),
    beta = as.numeric(beta_nfr),
    se = as.numeric(se_nfr),
    pval = as.numeric(p_nfr),
    effect_allele = as.character(ea),
    other_allele = as.character(oa),
    eaf = as.numeric(eaf),
    samplesize = as.numeric(n_nfr),
    id.exposure = "nfr",
    exposure = "nfr",
    stringsAsFactors = FALSE
  )

  out_df <- data.frame(
    SNP = as.character(snp),
    beta = as.numeric(beta_nuc),
    se = as.numeric(se_nuc),
    pval = as.numeric(p_nuc),
    effect_allele = as.character(ea),
    other_allele = as.character(oa),
    eaf = as.numeric(eaf),
    samplesize = as.numeric(n_nuc),
    id.outcome = "nuc",
    outcome = "nuc",
    stringsAsFactors = FALSE
  )

  out <- tryCatch({
    exp_dat <- format_data(
      exp_df,
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
      out_df,
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

    harm <- harmonise_data(exp_dat, out_dat, action = 2)
    if (nrow(harm) == 0) {
      return(list(steiger_dir = NA, steiger_p = NA_real_, r2_nfr = NA_real_, r2_nuc = NA_real_, note = "no_harmonized_row"))
    }

    dir <- directionality_test(harm)
    if (nrow(dir) == 0) {
      return(list(steiger_dir = NA, steiger_p = NA_real_, r2_nfr = NA_real_, r2_nuc = NA_real_, note = "no_directionality_row"))
    }

    dir_flag <- pick_first_col(dir, c("correct_causal_direction", "steiger_dir"))
    pval <- suppressWarnings(as.numeric(pick_first_col(dir, c("steiger_pval", "steiger.pval", "pval"))))
    r2_exp <- suppressWarnings(as.numeric(pick_first_col(dir, c("r2.exposure", "r2_exposure", "snp_r2.exposure", "rsq.exposure"))))
    r2_out <- suppressWarnings(as.numeric(pick_first_col(dir, c("r2.outcome", "r2_outcome", "snp_r2.outcome", "rsq.outcome"))))

    if (is.na(r2_exp)) {
      r_exp <- suppressWarnings(as.numeric(pick_first_col(dir, c("r.exposure", "r_exposure"))))
      if (!is.na(r_exp)) r2_exp <- r_exp ^ 2
    }
    if (is.na(r2_out)) {
      r_out <- suppressWarnings(as.numeric(pick_first_col(dir, c("r.outcome", "r_outcome"))))
      if (!is.na(r_out)) r2_out <- r_out ^ 2
    }

    list(
      steiger_dir = if (is.na(dir_flag)) NA else as.logical(dir_flag),
      steiger_p = pval,
      r2_nfr = r2_exp,
      r2_nuc = r2_out,
      note = "ok"
    )
  }, error = function(e) {
    list(steiger_dir = NA, steiger_p = NA_real_, r2_nfr = NA_real_, r2_nuc = NA_real_, note = paste0("error:", conditionMessage(e)))
  })

  out
}

rows <- vector("list", nrow(df))
for (i in seq_len(nrow(df))) {
  r <- df[i, ]

  t1 <- run_directionality_one(
    snp = r$hit1,
    ea = r$h1_effect_allele,
    oa = r$h1_other_allele,
    eaf = r$h1_eaf,
    beta_nfr = r$h1_nfr_beta,
    se_nfr = r$h1_nfr_se,
    p_nfr = r$h1_nfr_pval,
    n_nfr = r$samplesize_nfr,
    beta_nuc = r$h1_nuc_beta,
    se_nuc = r$h1_nuc_se,
    p_nuc = r$h1_nuc_pval,
    n_nuc = r$samplesize_nuc
  )

  t2 <- run_directionality_one(
    snp = r$hit2,
    ea = r$h2_effect_allele,
    oa = r$h2_other_allele,
    eaf = r$h2_eaf,
    beta_nfr = r$h2_nfr_beta,
    se_nfr = r$h2_nfr_se,
    p_nfr = r$h2_nfr_pval,
    n_nfr = r$samplesize_nfr,
    beta_nuc = r$h2_nuc_beta,
    se_nuc = r$h2_nuc_se,
    p_nuc = r$h2_nuc_pval,
    n_nuc = r$samplesize_nuc
  )

  rows[[i]] <- data.frame(
    nfr_peak = r$NFR_peak,
    nuc_peak = r$Nucleosomal_peak,
    hit1_snp = r$hit1,
    r2_nfr_hit1 = t1$r2_nfr,
    r2_nuc_hit1 = t1$r2_nuc,
    steiger_dir_hit1 = t1$steiger_dir,
    steiger_p_hit1 = t1$steiger_p,
    note_hit1 = t1$note,
    hit2_snp = r$hit2,
    r2_nfr_hit2 = t2$r2_nfr,
    r2_nuc_hit2 = t2$r2_nuc,
    steiger_dir_hit2 = t2$steiger_dir,
    steiger_p_hit2 = t2$steiger_p,
    note_hit2 = t2$note,
    stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, rows)

p_all <- c(out$steiger_p_hit1, out$steiger_p_hit2)
q_all <- rep(NA_real_, length(p_all))
valid <- which(!is.na(p_all) & is.finite(p_all) & p_all >= 0 & p_all <= 1)
if (length(valid) > 0) {
  qv <- tryCatch(qvalue(p_all[valid])$qvalues, error = function(e) p.adjust(p_all[valid], method = "BH"))
  q_all[valid] <- qv
}

n <- nrow(out)
out$steiger_q_hit1 <- q_all[seq_len(n)]
out$steiger_q_hit2 <- q_all[n + seq_len(n)]

q_thr <- as.numeric(opts$`q-threshold`)
assign_direction <- function(r2_nfr, r2_nuc, qv, qthr) {
  if (is.na(r2_nfr) || is.na(r2_nuc) || is.na(qv)) return("undetermined")
  if (qv < qthr && r2_nfr > r2_nuc) return("nfr_to_nuc")
  if (qv < qthr && r2_nuc > r2_nfr) return("nuc_to_nfr")
  "undetermined"
}

out$direction_hit1 <- mapply(assign_direction, out$r2_nfr_hit1, out$r2_nuc_hit1, out$steiger_q_hit1, MoreArgs = list(qthr = q_thr))
out$direction_hit2 <- mapply(assign_direction, out$r2_nfr_hit2, out$r2_nuc_hit2, out$steiger_q_hit2, MoreArgs = list(qthr = q_thr))

write.table(out, file = opts$out, sep = "\t", quote = FALSE, row.names = FALSE)
