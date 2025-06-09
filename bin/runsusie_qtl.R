#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(vcfR)
library(coloc)
library(dplyr)
library(stringr)
library(susieR)
library(optparse)

option_list <- list(
    make_option(c("--chr"), type = "character", help = "[Required] Looking at which chromosome"),
    make_option(c("--method"), type = "character", help = "[Required] Checking QTLs from which method"),
    make_option(c("--min_corr"), type = "numeric", help = "[Required] min abs corr between any snp pairs"),
    make_option(c("--num_L"), type = "numeric", help = "[Required] Max number of LD vars"),
    make_option(c("--outdir"), type = "character", help = "[Required] output directory")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# -----------until functions-----------

options(scipen = 100, digits = 4)
format_LD <- function(filename) {
    # tidy the sample list
    # samples to keep
    # only keep 281 samples - FUSION
    sample_info =  read.csv("../../data/sample_info/sample_level_covariates_atac_with_frag.tsv", sep = "\t")
    sample_info <- sample_info %>% filter(coarse_cluster_name == 'fusion.Type_1')
    # current list of samples - 284
    sample_list <- sample_info$SNG.1ST
    # three samples to exclude
    sample_exclude <- c("12004", "22011", "32071")
    sample_list <- sample_list[sample_list %in% sample_exclude == FALSE]
    # read vcf
    vcf = read.vcfR(filename)
    # get unique ids
    myID <- getID(vcf)
    vcf2 <- vcf[!duplicated(myID, incomparables = NA), ]
    dosages = extract.gt(vcf2, element="DS")
    dat <- mutate_all(as.data.frame(dosages), function(x) as.numeric(x))
    # keep the df with 281 samples
    idx <- match(sample_list, colnames(dat))
    dat_281 <- dat[,idx]
    #print(paste("Before removing duplicates, the number of snps inluced:", nrow(dat_281)))
    #dat_281 <- unique(dat_281)
    #print(paste("After removing duplicates, the number of snps inluced:", nrow(dat_281)))
    # transpose dataframe
    dat_281 <- t(dat_281)
    # turn into a dataframe
    dat_281 <- as.data.frame(dat_281)
    return(dat_281)
}

# input variables
chr = opts$chr
# the method
method = opts$method
# min abs correlation
min_corr = opts$min_corr
# number of var -L
num_L = opts$num_L
# outputs dir
outdir = opts$outdir

# -----------process inputs-----------
annotation = strsplit(method, split = "_counts")[[1]][1]
## annotations
NFR_sigQTLs_path = paste0("../../data/annotations/", annotation, "_sigQTL.bed")
## qtl nominal pass results
#NFR_qtl_summary_path = "../results/results_06062023/HMMRATAC_NFR_counts_multiotsu_normalized_F50_G5/HMMRATAC_NFR_counts_multiotsu_normalized_permutations_F50_G5.chr3.txt"
NFR_qtl_summary_path = paste0("../../results/results_06062023/", method, "/",
                              strsplit(method, split = "normalized")[[1]][1], 
                              "normalized_permutations", 
                              strsplit(method, split = "normalized")[[1]][2], 
                              ".", chr, ".txt")
n95cset = c()
#------------read in-----------------
NFR_sigQTLs = read.table(NFR_sigQTLs_path, col.names = c("Chr", "Start", "End"))
# filter sig nfrQTLs to the correct chromosome
NFR_sigQTLs_chr = NFR_sigQTLs %>% filter(Chr == chr)
NFR_qtl_summary = read.table(NFR_qtl_summary_path, 
                        sep = ' ', 
                        col.names = c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var', 'var_id', 
                                      'var_chr', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'slope', 'slope_se', 'best_hit'))
print(paste0("Creating susie models for sig QTLs from: ", method))
# start looping
for (i in 1:nrow(NFR_sigQTLs_chr)){
    # get the peak info - test with i=905 (nfr)
    chr = NFR_sigQTLs_chr[i,][['Chr']]
    start = NFR_sigQTLs_chr[i,][['Start']]
    end = NFR_sigQTLs_chr[i,][['End']]
    # create the qtl scan pheid
    pheid = paste0(chr, '.', start-1, '.', end)
    print(paste("Current looking at this sigQTL:", pheid))
    # peak df to coloc df
    ## nfr - based on pheid
    nfr_qtl_summary_peak = NFR_qtl_summary %>% filter(phe_id == pheid)
    ## coloc df - nfr
    ### reformat to coloc stats
    coloc_df1 <- data.frame(matrix(ncol = 0, nrow = nrow(nfr_qtl_summary_peak)))
    ### beta col
    coloc_df1$beta = nfr_qtl_summary_peak$slope
    ### varbeta col
    coloc_df1$varbeta = (nfr_qtl_summary_peak$slope_se)^2
    ### snp col
    coloc_df1$snp = nfr_qtl_summary_peak$var_id
    ### position col - excluding chr info now
    coloc_df1$position = nfr_qtl_summary_peak$var_from
    ### pvalues
    coloc_df1$pvalues = nfr_qtl_summary_peak$nom_pval
    # Get the LD matrix
    print("-------------Assembling LD matrix-------------")
    LD_chr = chr
    LD_start = coloc_df1 %>% slice_min(position) %>% select(position) %>% unique() %>% as.numeric()
    LD_end = coloc_df1 %>% slice_max(position) %>% select(position) %>% unique() %>% as.numeric()
    LD_pos = paste0(LD_chr, ":", LD_start-1, "-", LD_end+1)
    LD_file = paste0(LD_chr, "-", LD_start-1, "-", LD_end+1)
    command = "tabix -h ../../data/sample_info/fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz "
    command = paste0(command, LD_pos)
    command = paste0(command, " | bgzip > ../../data/sample_info/temp_vcfs/", LD_file, ".vcf.gz")
    print(command)
    system(command, intern=TRUE)
    # get the sample by var genotype info dat
    LD_dat = format_LD(paste0("../../data/sample_info/temp_vcfs/", LD_file, ".vcf.gz"))
    # filter to keep the same snps at the QTL merge df
    LD_dat = LD_dat[, which((names(LD_dat) %in% coloc_df1$snp)==TRUE)]
    # calculate cor
    LD_dat_m <- cor(LD_dat)
    # get coloc_d1
    print("-------------Assembling colco datalist -------------")
    # change the coloc dat to list
    coloc_d1 <- list()
    # beta
    coloc_d1$beta <- coloc_df1$beta
    # varbeta
    coloc_d1$varbeta <- coloc_df1$varbeta
    # add N
    coloc_d1$N = 281
    # sdY
    coloc_d1$sdY = 1.0
    # type col - quant now
    coloc_d1$type = "quant"
    # add LD
    coloc_d1$LD = LD_dat_m[coloc_df1$snp, coloc_df1$snp]
    # snp
    coloc_d1$snp <- coloc_df1$snp
    # position
    coloc_d1$position <- coloc_df1$position
    # pvalues
    coloc_d1$pvalues <- coloc_df1$pvalues
    # save dat
    save(coloc_d1, file = paste0(outdir, annotation, '_sigQTL_susie_models/', chr, '/', pheid, '_susie_dat.Rda'))
    # run susie
    print("-------------Running susie_rss -------------")
    S1 = susieR::susie_rss(z = coloc_d1$beta/sqrt(coloc_d1$varbeta), 
                       R = coloc_d1$LD, 
                       n = coloc_d1$N,
                       L = num_L,
                       estimate_residual_variance = TRUE,
                       estimate_prior_variance=TRUE,
                       min_abs_corr=min_corr
                       )
    # saving key outputs
    print("-------------Saving key outputs -------------")
    ## modified from Arushi's code
    ## susie obj
    # save the whole model to rda file
    save(S1, file = paste0(outdir, annotation, '_sigQTL_susie_models/', chr, '/', pheid, '.susie.Rda'))
    pipv = susie_get_pip(S1)
    pip_df = data.frame(list("pip"=pipv, "snp"=coloc_d1$snp))
    write.table(pip_df, file=paste0(outdir, annotation, '_sigQTL_susie_models/', chr, '/', pheid, '.susie_pipdf.tsv'),
                sep='\t', quote=F)
    # is there 95% credible set?
    ## how many 95% csets?
    n95 = length(S1$sets$cs)
    n95cset = c(n95cset, n95)
    ## best set
    best = row.names(S1$sets$purity[which.max(S1$sets$purity$min.abs.corr),])
    ## size of the main set
    len_1_95 = length(S1$sets$cs[[best]])
    ## save the main 95% cset
    write.table(pip_df[S1$sets$cs[[best]],], 
                file=paste0(outdir, annotation, '_sigQTL_susie_models/', chr, '/', pheid, '.susie_cset95.tsv'), 
                sep='\t', quote=F)

    print(paste(n95, "95% credible sets found. Length of the first 95% credible set =", len_1_95))
    print("-------------DONE for this QTL -------------\n")
}

NFR_sigQTLs_chr$n95cset = n95cset
write.table(NFR_sigQTLs_chr, 
            file=paste0(outdir, annotation, '_sigQTL_susie_models/', chr, '/', chr, '_susie_cset95_summary.tsv'), 
            sep='\t', quote=F)

print("\n-------------DONE for this chr -------------")

#command_done = paste0("vim ", outdir, annotation, "_sigQTL_susie_models/", annotation, "_", chr, "_done.txt")
#system(command_done, intern=TRUE)
