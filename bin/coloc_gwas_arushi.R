#!/usr/bin/env Rscript

# libraries
options(stringsAsFactors=FALSE)
library(coloc)
library(dplyr)
library(stringr)
library(susieR)
library(optparse)

option_list <- list(
    make_option(c("--gwas_dat"), type = "character", help = "[Required] Looking at which gwas susie model?"),
    make_option(c("--indir_qtl"), type = "character", help = "[Required] Location of the QTL summary"),
    make_option(c("--outdir_coloc"), type = "character", help = "[Required] colco results output directory"),
    make_option(c("--filter_window"), type = "numeric", help = "[Required] the mid-point distance between the peaks"),
    make_option(c("--region"), type = "character", help = "[Required] which type of QTL? NFR vs. Nucleosoaml")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# input variables
gwas_dat = opts$gwas_dat
# input dir QTL
indir_qtl = opts$indir_qtl
# outputs dir
outdir_coloc = opts$outdir_coloc
# filtering window size
filter_window = opts$filter_window
# qtl type
region = opts$region

# process input - GWAS side
# get the chr info -since QTLs based on chr
chr = strsplit(strsplit(gwas_dat, split='gwas_')[[1]][2], split="-")[[1]][1]
# gwas peak info
gwas_peak = strsplit(strsplit(gwas_dat, split='gwas_')[[1]][2], split="[.]")[[1]][1]

# process input - QTLs side
# susie models
NFR_models = paste0(indir_qtl, chr, "/", chr, "_susie_cset95_summary.tsv")
#Nucleosomal_models = paste0(indir_qtl, trait, "_Nucleosomal_sigQTL_susie_models/", chr, "/", chr, "_susie_cset95_summary.tsv")
# read in
NFR_models_summary = read.csv(NFR_models, sep ="\t", header=TRUE)
#Nucleosomal_models_summary = read.csv(Nucleosomal_models, sep ="\t")
# if no sigqtls
if ((nrow(NFR_models_summary) == 0)){
    # track lists
    df_res_nfr = as.data.frame(matrix(nrow=0,ncol=12))
    options(scipen=999)
    save(df_res_nfr, file = paste0(outdir_coloc, "gwas_", gwas_peak, "_res_summary_", filter_window, "win.Rda"))
    print("NO significant qtls in this chr")
    quit()
}
# filter to those with >=1 n95cset
NFR_models_summary_sig = NFR_models_summary %>% filter(n95cset >= 1)
#Nucleosomal_models_summary_sig = Nucleosomal_models_summary %>% filter(n95cset >= 1)
# get a peak mid point column
NFR_models_summary_sig$Mid = 0.5*(NFR_models_summary_sig$Start+NFR_models_summary_sig$End)
#Nucleosomal_models_summary_sig$Mid = 0.5*(Nucleosomal_models_summary_sig$Start+Nucleosomal_models_summary_sig$End)

# gwas peak mid point
gwas_mid = strtoi(strsplit(gwas_peak, split = "-")[[1]][2])+500000

print(paste0("---------GWAS and ", region, "QTLs---------"))
# track lists - nfr
df_res_nfr = as.data.frame(matrix(nrow=0,ncol=12))
nearby_peaks = NFR_models_summary_sig %>% 
        filter(abs(Mid - gwas_mid) <= filter_window)
# loop through nearby peaks
if (nrow(nearby_peaks) > 0){
    for (i in 1:nrow(nearby_peaks)){
        print(paste0("---------Studying coloc between GWAS and ", region, "QTLs---------"))
        print(gwas_peak)
        # get nfr peak info
        nfr_start = nearby_peaks[i,][['Start']]
        nfr_end = nearby_peaks[i,][['End']]
        # build the nfr pheid
        nfr_pheid = paste0(chr, '.', nfr_start-1, '.', nfr_end)
        print(paste(region, nfr_pheid))
        S_nfr_path = paste0(indir_qtl, chr, "/", nfr_pheid, ".susie.rda")
        # load the susie models
        load(gwas_dat)
        S_gwas = S1
        load(S_nfr_path)
        S_nfr = S1
        # run coloc
        susie.res=coloc.susie(S_gwas, S_nfr)
        # save outputs
        ## res
        save(susie.res, 
            file = paste0(outdir_coloc, 
                        gwas_peak, "_", 
                        nfr_pheid, "_res.Rda"))
        ## pph4
        summary_df = susie.res$summary
        if(is.null(summary_df) == FALSE){
            summary_df$QTL_peak = rep(nfr_pheid, nrow(summary_df))
            summary_df$gwas_peak = rep(gwas_peak, nrow(summary_df))
            # bind
            df_res_nfr = rbind(df_res_nfr, summary_df)
        }
    }
}
# save summary
options(scipen=999)
save(df_res_nfr, file = paste0(outdir_coloc, "gwas_", gwas_peak, "_res_summary_", filter_window, "win.Rda"))


# if (region == "NFR"){
#     print("---------GWAS and nfrQTLs---------")
#     # track lists - nfr
#     df_res_nfr = as.data.frame(matrix(nrow=0,ncol=12))
#     nearby_peaks = NFR_models_summary_sig %>% 
#             filter(abs(Mid - gwas_mid) <= filter_window)
#     # loop through nearby peaks
#     if (nrow(nearby_peaks) > 1){
#         for (i in 1:nrow(nearby_peaks)){
#             print("---------Studying coloc between GWAS and nfrQTLs---------")
#             print(gwas_peak)
#             # get nfr peak info
#             nfr_start = nearby_peaks[i,][['Start']]
#             nfr_end = nearby_peaks[i,][['End']]
#             # build the nfr pheid
#             nfr_pheid = paste0(chr, '.', nfr_start-1, '.', nfr_end)
#             print(paste(trait, nfr_pheid))
#             S_nfr_path = paste0(indir_qtl, trait, "_NFR_sigQTL_susie_models/", chr, "/", nfr_pheid, ".susie.Rda")
#             # load the susie models
#             load(gwas_dat)
#             S_gwas = S1
#             load(S_nfr_path)
#             S_nfr = S1
#             # run coloc
#             susie.res=coloc.susie(S_gwas, S_nfr)
#             # save outputs
#             ## res
#             save(susie.res, 
#                 file = paste0(outdir_coloc, "coloc_", trait, "_NFR", "/", GWAS_annot, "/", 
#                             gwas_peak, "_", 
#                             nfr_pheid, "_res.Rda"))
#             ## pph4
#             summary_df = susie.res$summary
#             if(is.null(summary_df) == FALSE){
#                 summary_df$NFR_peak = rep(nfr_pheid, nrow(summary_df))
#                 summary_df$gwas_peak = rep(gwas_peak, nrow(summary_df))
#                 # bind
#                 df_res_nfr = rbind(df_res_nfr, summary_df)
#             }
#         }
#     }
#     # save summary
#     options(scipen=999)
#     save(df_res_nfr, file = paste0(outdir_coloc, "coloc_", trait, "_NFR", "/", GWAS_annot, "/", "summary/gwas_", gwas_peak, "_res_summary_", filter_window, "win.Rda"))
# }

# if (region == "Nucleosomal"){
#     print("---------GWAS and nucQTLs---------")
#     # track lists - nuc
#     df_res_nuc = as.data.frame(matrix(nrow=0,ncol=12))
#     nearby_peaks = Nucleosomal_models_summary_sig %>% 
#             filter(abs(Mid - gwas_mid) <= filter_window)
#     # loop through nearby peaks
#     if (nrow(nearby_peaks) > 1){
#         for (i in 1:nrow(nearby_peaks)){
#             print("---------Studying coloc between GWAS and nucQTLs---------")
#             print(gwas_peak)
#             # get nuc peak info
#             nuc_start = nearby_peaks[i,][['Start']]
#             nuc_end = nearby_peaks[i,][['End']]
#             # build the nuc pheid
#             nuc_pheid = paste0(chr, '.', nuc_start-1, '.', nuc_end)
#             print(paste(trait, nuc_pheid))
#             S_nuc_path = paste0(indir_qtl, trait, "_Nucleosomal_sigQTL_susie_models/", chr, "/", nuc_pheid, ".susie.Rda")
#             # load the susie models
#             load(gwas_dat)
#             S_gwas = S1
#             load(S_nuc_path)
#             S_nuc = S1
#             # run coloc
#             susie.res=coloc.susie(S_gwas, S_nuc)
#             # save outputs
#             ## res
#             save(susie.res, 
#                 file = paste0(outdir_coloc, "coloc_", trait, "_Nucleosomal", "/", GWAS_annot, "/", 
#                             gwas_peak, "_", 
#                             nuc_pheid, "_res.Rda"))
#             ## pph4
#             summary_df = susie.res$summary
#             if(is.null(summary_df) == FALSE){
#                 summary_df$Nucleosomal_peak = rep(nuc_pheid, nrow(summary_df))
#                 summary_df$gwas_peak = rep(gwas_peak, nrow(summary_df))
#                 # bind
#                 df_res_nuc = rbind(df_res_nuc, summary_df)
#             }
#         }
#     }

#     # write table
#     # save summary
#     options(scipen=999)
#     save(df_res_nuc, file = paste0(outdir_coloc, "coloc_", trait, "_Nucleosomal", "/", GWAS_annot, "/", "summary/gwas_", gwas_peak, "_res_summary_", filter_window, "win.Rda"))
# }
  