#!/usr/bin/env Rscript

# libraries
options(stringsAsFactors=FALSE)
library(coloc)
library(dplyr)
library(stringr)
library(susieR)
library(optparse)

option_list <- list(
    make_option(c("--chr"), type = "character", help = "[Required] Looking at which chromosome"),
    make_option(c("--nfr"), type = "character", help = "[Required] Nfr susie summary dir"),
    make_option(c("--nuc"), type = "character", help = "[Required] Nuc susie summary dir"),
    make_option(c("--indir_nfr"), type = "character", help = "[Required] Location of the nfr susie models"),
    make_option(c("--indir_nuc"), type = "character", help = "[Required] Location of the nuc susie models"),
    make_option(c("--outdir_coloc"), type = "character", help = "[Required] colco results output directory"),
    make_option(c("--filter_window"), type = "numeric", help = "[Required] the mid-point distance between the peaks")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# input variables
chr = opts$chr
# input dir
indir_nfr = opts$indir_nfr
indir_nuc = opts$indir_nuc
# outputs dir
outdir_coloc = opts$outdir_coloc
# filtering window size
filter_window = opts$filter_window

# process input
# susie models
#NFR_models = paste0(indir, trait, "_NFR_sigQTL_susie_models/", chr, "/", chr, "_susie_cset95_summary.tsv")
#Nucleosomal_models = paste0(indir, trait, "_Nucleosomal_sigQTL_susie_models/", chr, "/", chr, "_susie_cset95_summary.tsv")
NFR_models = opts$nfr
Nucleosomal_models = opts$nuc
# read in
NFR_models_summary = read.csv(NFR_models, sep ="\t", header = TRUE)
Nucleosomal_models_summary = read.csv(Nucleosomal_models, sep ="\t", header = TRUE)
# if no sigqtls
if ((nrow(NFR_models_summary) == 0) || (nrow(Nucleosomal_models_summary) == 0) ){
    print("NO significant qtls in this chr")
    # track lists
    df_res = as.data.frame(matrix(nrow=0,ncol=12))
    save(df_res, file = paste0(outdir_coloc, chr, "_res_summary_", filter_window, "win.Rda"))
    quit()
}
# filter to those with >=1 n95cset
NFR_models_summary_sig = NFR_models_summary %>% filter(n95cset >= 1)
Nucleosomal_models_summary_sig = Nucleosomal_models_summary %>% filter(n95cset >= 1)
# get a peak mid point column
NFR_models_summary_sig$Mid = 0.5*(NFR_models_summary_sig$Start+NFR_models_summary_sig$End)
Nucleosomal_models_summary_sig$Mid = 0.5*(Nucleosomal_models_summary_sig$Start+Nucleosomal_models_summary_sig$End)

# track lists
df_res = as.data.frame(matrix(nrow=0,ncol=12))
# loop through the peaks
for (i in 1:nrow(NFR_models_summary_sig)){
#for (i in 1:10){
    # get the NFR peak info
    nfr_start = NFR_models_summary_sig[i,][['Start']]
    nfr_end = NFR_models_summary_sig[i,][['End']]
    # get the peak mid point location
    nfr_mid = NFR_models_summary_sig[i,][['Mid']]
    # create the qtl scan pheid
    nfr_pheid = paste0(chr, '.', nfr_start-1, '.', nfr_end)
    # check if there is nearby nucleo peak
    nearby_peaks = Nucleosomal_models_summary_sig %>% 
        filter(abs(Mid - nfr_mid) <= filter_window)
    # if it is not empty
    if (nrow(nearby_peaks) > 0){
        # loop through the nearby peaks df
        for (j in 1:nrow(nearby_peaks)){
            # get nucleo peak info
            nucleosomal_start = nearby_peaks[j,][['Start']]
            nucleosomal_end = nearby_peaks[j,][['End']]
            # build the nucleosomal pheid
            nucleosomal_pheid = paste0(chr, '.', nucleosomal_start-1, '.', nucleosomal_end)
            print("---------Studying coloc between those two pairs---------")
            print(paste("NFR: ", nfr_pheid))
            print(paste("Nuleosomal: ", nucleosomal_pheid))
            # get the susie models
            S_nfr_path = paste0(indir_nfr, nfr_pheid, ".susie.rda")
            S_nucleosomal_path = paste0(indir_nuc, nucleosomal_pheid, ".susie.rda")
            # load the models
            load(file = S_nfr_path)
            print(S_nfr_path)
            S_nfr = S1
            load(file = S_nucleosomal_path)
            print(S_nucleosomal_path)
            S_nucleosomal = S1
            # run coloc
            susie.res=coloc.susie(S_nfr,S_nucleosomal)
            # save outputs
            ## res
            save(susie.res, file = paste0(outdir_coloc, nfr_pheid, "_", nucleosomal_pheid, "_res.Rda"))
            ## pph4
            summary_df = susie.res$summary
            if(is.null(summary_df) == FALSE){
                summary_df$NFR_peak = rep(nfr_pheid, nrow(summary_df))
                summary_df$Nucleosomal_peak = rep(nucleosomal_pheid, nrow(summary_df))
                # bind
                df_res = rbind(df_res, summary_df)
            }
            
        }
        #break

    }
    
    #break
}
# write table
# save output
save(df_res, file = paste0(outdir_coloc, chr, "_res_summary_", filter_window, "win.Rda"))