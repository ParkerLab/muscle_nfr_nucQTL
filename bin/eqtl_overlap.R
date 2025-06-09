#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)
library(dplyr)
library(stringr)
library(susieR)

option_list <- list(
    make_option(c("--celltype"), type = "character", help = "[Required] The celltype that we are checking"),
    make_option(c("--eqtl"), type = "character", help = "[Required] eQTL summary saved location"),
    make_option(c("--output"), type = "character", help = "[Required] eQTL overlap info saved output filename")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

celltype = opts$celltype

# function
nfr_e_overlap <- function(celltype, Chr, eqtl_summary_df, window){
    # read in
    # process input
    nfrqtl_susie_summary_path = paste0('../../results/results_07162024/', celltype, '_nfr_susie/', Chr, '/', Chr, '_susie_cset95_summary.tsv')
    nucqtl_susie_summary_path = paste0('../../results/results_07162024/', celltype, '_nuc_susie/', Chr, '/', Chr, '_susie_cset95_summary.tsv')
    nfrqtl_summary_df = read.table(nfrqtl_susie_summary_path)
    nucqtl_summary_df = read.table(nucqtl_susie_summary_path)
    overlap_list = c()
    for (i in 1:nrow(nfrqtl_summary_df)){
        # process info
        Chr = nfrqtl_summary_df[i,]$Chr
        start = max(0, nfrqtl_summary_df[i,]$Start-window)
        end = nfrqtl_summary_df[i,]$End+window
        #print(nfrqtl_summary_df[i,]$n95cset)
        if (nfrqtl_summary_df[i,]$n95cset == 0){
            overlap_list = c(overlap_list, 0)
        } else {
            if(nfrqtl_summary_df[i,]$n95cset == 1){
                # read in the nfrqtl_cset
                nfrcset_path = paste0('../../results/results_07162024/', celltype, '_nfr_susie/', Chr, '/', 
                                      Chr, '.', nfrqtl_summary_df[i,]$Start-1, '.', 
                                      nfrqtl_summary_df[i,]$End, '.susie_cset95.tsv')
                nfr_cset = read.table(nfrcset_path)
                # get the nfr snps
                nfr_snps = nfr_cset$snp
            } else {
                #print(nfrqtl_summary_df[i,])
                #print(nfrqtl_summary_df[i,]$n95cset)
                # load the susie result
                susie_path = paste0('../../results/results_07162024/', celltype, '_nfr_susie/', Chr, '/', 
                                    Chr, '.', nfrqtl_summary_df[i,]$Start-1, '.', 
                                    nfrqtl_summary_df[i,]$End, '.susie.rda')
                load(susie_path)
                # get pipdf
                pip_path = paste0('../../results/results_07162024/', celltype, '_nfr_susie/', Chr, '/', 
                                    Chr, '.', nfrqtl_summary_df[i,]$Start-1, '.', 
                                    nfrqtl_summary_df[i,]$End, '.susie_pipdf.tsv')
                pip_df = read.table(pip_path)
                snp_list = c()
                for(i in 1:length(S1$sets$cs)){
                    snp_list = c(snp_list, S1$sets$cs[[i]])
                }
                # get the snps
                snp_df = pip_df[snp_list,]
                # get nfr snps
                nfr_snps = snp_df$snp
                # save this df
                fullsnp_path = paste0('../../results/results_07162024/', celltype, '_nfr_susie/', Chr, '/', 
                                      Chr, '.', nfrqtl_summary_df[i,]$Start-1, '.', 
                                      nfrqtl_summary_df[i,]$End, '.susie_cset95_full.tsv')
                write.table(snp_df, file=fullsnp_path, 
                            sep='\t', quote=F)
            }
            # filter for possible hits in eqtl
            eqtl_tmp_df = eqtl_summary_df %>% 
                            filter(chr == Chr) %>% 
                            filter(pos >= start) %>%
                            filter(pos <= end)
            snp_list = c()
            if (nrow(eqtl_tmp_df) == 0){
                overlap_list = c(overlap_list, 0)
            } else {
                for (j in 1:nrow(eqtl_tmp_df)){
                    snp_df = read.table(eqtl_tmp_df[j,]$path)
                    snp_list = c(snp_list, snp_df$V7)
                }
                # check if overlap
                if (length(intersect(unlist(nfr_snps), unlist(snp_list))) > 0){
                    overlap_list = c(overlap_list, 1)
                } else {
                    overlap_list = c(overlap_list, 0)
                }
            }
        
        } 

    }
    nfrqtl_summary_df$eqtl_overlap = overlap_list
    
    return(nfrqtl_summary_df)
}


eqtl_summary_df = read.table(opts$eqtl)
# init empty df
nfrqtl_summary_df_all = data.frame()
for (i in 1:22){
    Chr = paste0('chr', i)
    print(Chr)
    # apply function
    nfrqtl_summary_df = nfr_e_overlap(celltype, Chr, eqtl_summary_df, 250000)
    # nuc coloc info
    summary_path = paste0('../../results/results_07172024/coloc_', celltype, 
                            '/chr', i, '/chr', i, '_res_summary_5000win.Rda')
    # get nfrpeak name
    nfrqtl_summary_df$NFR_peak = paste0(nfrqtl_summary_df$Chr, '.',
                                    nfrqtl_summary_df$Start-1, '.',
                                    nfrqtl_summary_df$End)
    # load
    load(summary_path)
    # if empty
    if (nrow(df_res) > 0) {
        # filter for pph4>0.8
        df_coloc = df_res %>% filter(PP.H4.abf > 0.8)
        # assign 0/1
        df_coloc$nuc_coloc=1
        df_coloc = df_coloc %>% select(c('NFR_peak', 'nuc_coloc'))

        # merge
        nfrqtl_summary_df <- nfrqtl_summary_df %>% merge(df_coloc, by='NFR_peak', all.x=T)
        # asign 0
        nfrqtl_summary_df[is.na(nfrqtl_summary_df)] <- 0
        
    } else {
        nfrqtl_summary_df$nuc_coloc=0
    }
    # rbind
    nfrqtl_summary_df_all = rbind(nfrqtl_summary_df_all, nfrqtl_summary_df)
}

# save results
filename = opts$output
write.table(nfrqtl_summary_df_all, file=filename, sep = '\t', quote=F)