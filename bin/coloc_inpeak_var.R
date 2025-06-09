#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)
library(coloc)
library(dplyr)
library(stringr)
library(susieR)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(patchwork)
library(Rfast)

option_list <- list(
    make_option(c("--chr"), type = "character", help = "[Required] Chromosome"),
    make_option(c("--celltype"), type = "character", help = "[Required] Celltype"),
    make_option(c("--gwas_summary"), type = "character", help = "[Required] Summary file per gwas signal coloc results"),
    make_option(c("--nominal_path"), type = "character", help = "[Required] Where is the directory for nfrqtl nominal results"),
    make_option(c("--nfrqtl_susie_dir"), type = "character", help = "[Required] Directory of the NFR QTL susie models"),
    make_option(c("--gwas_susie_dir"), type = "character", help = "[Required] Directory of the GWAS susie models"),
    make_option(c("--output"), type = "character", help = "[Required] Summary file for coloc gwas inpeak information")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)



nominal_path = opts$nominal_path
nfrqtl_susie_dir = opts$nfrqtl_susie_dir
gwas_susie_dir = opts$gwas_susie_dir
gwas_summary = opts$gwas_summary
output = opts$output
chr = opts$chr
celltype = opts$celltype

df = read.table(gwas_summary, sep = ' ')
# use high PPH4 threshold
df = df %>% filter(PP.H4.abf > 0.9)
# add chr col
df$Chr <- apply(df, 1, function(row) strsplit(row['gwas_peak'], '-')[[1]][1])
# filter for chr
tmp_df <- df %>% filter(Chr == chr) %>% filter(Celltype == celltype)
# check if empty
if(nrow(tmp_df) == 0){
    print("No GWAS signals for this chr and celltype")
    inpeak_df <- data.frame(
                            file_qtls = character(),
                            file_gwass = character(),
                            hit_vars = character(),
                            pph4s = numeric(),
                            stringsAsFactors = FALSE  # Keep character vectors as is, do not convert to factors
                            )

    write.table(inpeak_df, quote = F, sep = '\t', file=output)
    quit()
}

# read in nominal pass results
# nominal_pass result
path_nfr <- paste0(nominal_path, celltype, '_nfr/nominal_', chr, '.txt')
# read in nominal result
df_nfr <- read.table(path_nfr, sep = ' ',
                    col.names = c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 
                                    'n_var_in_cis', 'dist_phe_var', 'var_id',
                                    'var_chr', 'var_from', 'var_to', 'nom_pval', 
                                    'r_squared', 'slope', 'slope_se', 'best_hit'))

file_qtls = c()
file_gwass = c()
hit_vars = c()
sources = c()
pph4s = c()

# Assume df is your data frame already loaded in R
for (i in 1:nrow(tmp_df)) {
    print(i)
    # get gwas_peak info
    gwas_peak <- tmp_df$gwas_peak[i]
    print(gwas_peak)
    # qtl peak info
    qtl_peak <- tmp_df$QTL_peak[i]
    print(qtl_peak)
    # get chr information
    #chr <- strsplit(gwas_peak, '-')[[1]][1]

    # get celltype information
    #celltype <- df$Celltype[i]

    # get credible set information
    #gwas_cs95_df <- read.table(paste0(gwas_susie_dir, 'gwas_', gwas_peak, '.susie_cset95.tsv'))
    #qtl_cs95_df <- read.table(paste0(nfrqtl_susie_dir, celltype, '_nfr_susie/', chr, '/', qtl_peak, '.susie_cset95.tsv'))

    # count either gwas or qtl cs snps as possible source
    possible_vars = c(tmp_df[i,]$hit1, tmp_df[i,]$hit2) %>% unique()

    target_df <- df_nfr %>%
        filter(phe_id == qtl_peak & var_id %in% possible_vars) %>%
        filter(dist_phe_var==0)
    #print(nrow(target_df))
    if(nrow(target_df) == 0){
        next
    }
    for(j in nrow(target_df)){
        file_qtl = paste0(nfrqtl_susie_dir, celltype, '_nfr_susie/', chr, '/', qtl_peak, '_susie_dat.Rda')
        file_gwas = paste0(gwas_susie_dir, 'gwas_', gwas_peak, '_susie_dat.Rda')
        # get hit var
        hit_var = target_df[j,]$var_id
        # pph4
        pph4 = tmp_df[i,]$PP.H4.abf
        # append
        file_qtls = c(file_qtls, file_qtl)
        file_gwass = c(file_gwass, file_gwas)
        hit_vars = c(hit_vars, hit_var)
        print(hit_var, length(hit_vars))
        pph4s = c(pph4s, pph4)
    }
}
# Create a data frame
inpeak_df <- data.frame(
  file_qtls = file_qtls,
  file_gwass = file_gwass,
  hit_vars = hit_vars,
  pph4s = pph4s,
  stringsAsFactors = FALSE  # Keep character vectors as is, do not convert to factors
)
# Create a data frame
inpeak_df <- data.frame(
  file_qtls = file_qtls,
  file_gwass = file_gwass,
  hit_vars = hit_vars,
  pph4s = pph4s,
  stringsAsFactors = FALSE  # Keep character vectors as is, do not convert to factors
)

write.table(inpeak_df, quote = F, sep = '\t', file=output)