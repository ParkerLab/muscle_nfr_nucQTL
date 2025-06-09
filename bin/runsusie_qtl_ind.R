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
    make_option(c("--bed"), type = "character", help = "[Required] directory of corrected bed file"),
    make_option(c("--qtl"), type = "character", help = "[Required] directory of fdr corrected qtl results"),
    make_option(c("--nominal"), type = "character", help = "[Required] directory of nominal pass qtl results"),
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
    sample_list <- sample_info$SAMPLE
    # three samples to exclude - deindentified
    sample_exclude <- c("donor1", "donor2", "donor3")
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
# min abs correlation
min_corr = opts$min_corr
# number of var -L
num_L = opts$num_L
# outputs dir
outdir = opts$outdir
# dir for corrected bed files
bed = opts$bed
# fdr corrected qtl results
qtl = opts$qtl
# nominal pass result
nominal = opts$nominal

# -----------process inputs-----------
# check and read
## count mat locations
dat = as.data.frame(read.table(bed, comment.char=''), header = TRUE, sep="\t",stringsAsFactors=FALSE, quote="")
# assgin headers
names(dat) <- dat[1,]
dat <- dat[-1,]
count_mat = dat                     
n95cset = c()
#------------read in-----------------
fdr_qtls = read.csv(qtl)
# filter sig nfrQTLs to the correct chromosome
NFR_sigQTLs_chr = fdr_qtls %>% 
    filter(phe_chr == chr) %>% 
    filter(adj_beta_qval < 0.05) %>% 
    select(c('phe_chr', 'phe_from', 'phe_to')) %>%
    rename(Chr=phe_chr, Start=phe_from, End=phe_to)

NFR_qtl_summary = read.table(nominal, 
                        sep = ' ', 
                        col.names = c('phe_id', 'phe_chr', 'phe_from', 'phe_to', 'phe_strd', 'n_var_in_cis', 'dist_phe_var', 'var_id', 
                                      'var_chr', 'var_from', 'var_to', 'nom_pval', 'r_squared', 'slope', 'slope_se', 'best_hit'))
print(paste0("Creating susie models for sig QTLs from: ", qtl))
# start looping
if (nrow(NFR_sigQTLs_chr) != 0){
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
    #if (nrow(nfr_qtl_summary_peak) == 0) {
    #    n95cset = c(n95cset, 0)
    #    next
    #}
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
    # get y
    y = count_mat %>% 
        filter(id == pheid) %>% 
        select(-c('#chr', 'start', 'end', 'id', 'gid', 'strd')) %>% 
        t() %>% 
        c() %>% 
        as.numeric()
    # get coloc_d1
    print("-------------Assembling coloc datalist -------------")
    # filter colocdf1
    coloc_df1 = coloc_df1 %>% filter(snp %in% rownames(LD_dat_m))
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
    coloc_d1$snp <- rownames(LD_dat_m)
    # position
    coloc_d1$position <- coloc_df1$position
    # pvalues
    coloc_d1$pvalues <- coloc_df1$pvalues
    # also y
    coloc_d1$y <- y
    # save dat
    save(coloc_d1, file = paste0(outdir, pheid, '_susie_dat.Rda'))
    # run susie
    print("-------------Running susie using individual level data -------------")
    S1 = susie(data.matrix(LD_dat), y, L = num_L,
               estimate_residual_variance = TRUE, 
               estimate_prior_variance = TRUE,
               verbose = TRUE,
               min_abs_cor = min_corr)
    # saving key outputs
    print("-------------Saving key outputs -------------")
    ## modified from Arushi's code
    ## susie obj
    # save the whole model to rda file
    save(S1, file = paste0(outdir, pheid, '.susie.rda'))
    pipv = susie_get_pip(S1)
    pip_df = data.frame(list("pip"=pipv, "snp"=coloc_d1$snp))
    write.table(pip_df, file=paste0(outdir, pheid, '.susie_pipdf.tsv'),
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
                file=paste0(outdir, pheid, '.susie_cset95.tsv'), 
                sep='\t', quote=F)

    print(paste(n95, "95% credible sets found. Length of the first 95% credible set =", len_1_95))
    print("-------------DONE for this QTL -------------")
}
NFR_sigQTLs_chr$n95cset = n95cset
write.table(NFR_sigQTLs_chr, 
            file=paste0(outdir, chr, '_susie_cset95_summary.tsv'), 
            sep='\t', quote=F)

print("\n-------------DONE for this chr -------------")
}
# if no sigqtl in this chr
if (nrow(NFR_sigQTLs_chr) == 0) {
    write.table(NFR_sigQTLs_chr, 
            file=paste0(outdir, chr, '_susie_cset95_summary.tsv'), 
            sep='\t', quote=F)

    print("\n-------------no sigqtl in this chr -------------")
}



#command_done = paste0("vim ", outdir, annotation, "_sigQTL_susie_models/", annotation, "_", chr, "_done.txt")
#system(command_done, intern=TRUE)
