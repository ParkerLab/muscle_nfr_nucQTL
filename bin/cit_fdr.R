#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)
library(cit)
library(dplyr)
library(vcfR)
library(stringr)

option_list <- list(
    make_option(c("--celltype"), type = "character", help = "[Required] Which celltype?"),
    make_option(c("--cit_outdir"), type = "character", help = "[Required] Cit permutations results saved location"),
    make_option(c("--PPH4_threshold"), type = "numeric", help = "[Required] Filtering threshold for PPh4 value"),
    make_option(c("--fdr_outdir"), type = "character", help = "[Required] FDR corrected cir results saved location")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

celltype = opts$celltype
cit_outdir = opts$cit_outdir
PPH4_threshold = opts$PPH4_threshold
fdr_outdir = opts$fdr_outdir
#print(trait)
cit_results_path = paste0(cit_outdir, "cit_PPH4_threshold", PPH4_threshold, "/")
# process inputs
l = list.files(path = cit_results_path)
# load the perm results - nocov
idx = str_detect(l, paste0("^", celltype))
trait_results_nocov = l[unlist(idx)]
# init trait results
trait_H1_results = list()
trait_H2_results = list()
trait_H3_results = list()
trait_H4_results = list()
# init pair info lists
hit1 = c()
hit2 = c()
nfr_peak = c()
nuc_peak = c()
for(i in 1:length(trait_results_nocov)){
    cit_result_file = paste0(cit_results_path, trait_results_nocov[i])
    #print(cit_result_file)
    # load the results
    load(cit_result_file)
    trait_H1_results[[i]] = results$H1
    trait_H2_results[[i]] = results$H2
    trait_H3_results[[i]] = results$H3
    trait_H4_results[[i]] = results$H4
    # get the pair info
    hit1 = c(hit1, strsplit(trait_results_nocov[i], split= "-")[[1]][2])
    hit2 = c(hit2, strsplit(trait_results_nocov[i], split= "-")[[1]][3])
    nfr_peak = c(nfr_peak, strsplit(trait_results_nocov[i], split= "-")[[1]][4])
    nuc_peak = c(nuc_peak, strsplit(trait_results_nocov[i], split= "-")[[1]][5])
}
# get the cit pair info
pair_info_df = data.frame(hit1 = hit1, 
                          hit2 = hit2,
                          nfr_peak = nfr_peak,
                          nuc_peak = nuc_peak)
# run fdr
## H1
fdr_results_H1 = fdr.cit(trait_H1_results, cl=.95, c1=NA)
fdr_info_df_H1 = cbind(pair_info_df, fdr_results_H1)
## H2
fdr_results_H2 = fdr.cit(trait_H2_results, cl=.95, c1=NA)
fdr_info_df_H2 = cbind(pair_info_df, fdr_results_H2)
## H3
fdr_results_H3 = fdr.cit(trait_H3_results, cl=.95, c1=NA)
fdr_info_df_H3 = cbind(pair_info_df, fdr_results_H3)
## H4
fdr_results_H4 = fdr.cit(trait_H4_results, cl=.95, c1=NA)
fdr_info_df_H4 = cbind(pair_info_df, fdr_results_H4)

# save results
write.table(fdr_info_df_H1, 
            file=paste0(fdr_outdir, celltype, '_citfdr_cov_', PPH4_threshold, '_H1_summary.tsv'), 
            sep='\t', quote=F)
write.table(fdr_info_df_H2, 
            file=paste0(fdr_outdir, celltype, '_citfdr_cov_', PPH4_threshold, '_H2_summary.tsv'), 
            sep='\t', quote=F)
write.table(fdr_info_df_H3, 
            file=paste0(fdr_outdir, celltype, '_citfdr_cov_', PPH4_threshold, '_H3_summary.tsv'), 
            sep='\t', quote=F)
write.table(fdr_info_df_H4, 
            file=paste0(fdr_outdir, celltype, '_citfdr_cov_', PPH4_threshold, '_H4_summary.tsv'), 
            sep='\t', quote=F)