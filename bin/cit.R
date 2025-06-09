#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(optparse)
library(cit)
library(dplyr)
library(vcfR)

option_list <- list(
    make_option(c("--pair"), type = "character", help = "[Required] Cit pair info order: trait_hit1_hit2_nfrpeak_nucpeak"),
    make_option(c("--outdir"), type = "character", help = "[Required] Cit permutations results saved location"),
    make_option(c("--PPH4_threshold"), type = "numeric", help = "[Required] Filtering threshold for PPh4 value"),
    make_option(c("--cov"), type = "character", help = "[Required] location of the technical + genoPC cov file location"),
    make_option(c("--count_mat_dir"), type = "character", help = "[Required] location of the count/score matrix")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

cit_pair = opts$pair
outdir = opts$outdir
PPH4_threshold = opts$PPH4_threshold
cov_file = opts$cov

# process information
celltype = strsplit(cit_pair, split = "-")[[1]][1]
hit1 = strsplit(cit_pair, split = "-")[[1]][2]
hit2 = strsplit(cit_pair, split = "-")[[1]][3]
nfr_peak = strsplit(cit_pair, split = "-")[[1]][4]
nuc_peak = strsplit(cit_pair, split = "-")[[1]][5]
print(cit_pair)
print(paste(celltype, hit1, hit2, nfr_peak, nuc_peak))
print(nuc_peak)
# get sub_vcf
## chr info
chr = strsplit(nfr_peak, split = "[.]")[[1]][1]
## vcf start position
vcf_start_pos = min(strsplit(nfr_peak, split = "[.]")[[1]][2], 
                    strsplit(nfr_peak, split = "[.]")[[1]][3],
                    strsplit(nuc_peak, split = "[.]")[[1]][2], 
                    strsplit(nuc_peak, split = "[.]")[[1]][3])
vcf_start_pos = strtoi(vcf_start_pos) - 500000
## vcf end position
vcf_end_pos = max(strsplit(nfr_peak, split = "[.]")[[1]][2], 
                  strsplit(nfr_peak, split = "[.]")[[1]][3],
                  strsplit(nuc_peak, split = "[.]")[[1]][2], 
                  strsplit(nuc_peak, split = "[.]")[[1]][3])
vcf_end_pos = strtoi(vcf_end_pos) + 500000
## vcf slicing info
vcf_pos = paste0(chr, ":", vcf_start_pos, "-", vcf_end_pos)
command = "tabix -h ../../data/sample_info/fusion.filtered-vcf.maf0.05-hwe1e6.vcf.gz "
command = paste0(command, vcf_pos)
command = paste0(command, " | bgzip > ../../data/sample_info/temp_vcfs/", vcf_pos, ".vcf.gz")
print(command)
system(command, intern=TRUE)

# --------------------SNP info-----------------------
# get the vcf file
filename = paste0("../../data/sample_info/temp_vcfs/", vcf_pos, ".vcf.gz")
# tidy the dosage information
# tidy the sample list
# samples to keep
# only keep 281 samples - FUSION
sample_info =  read.csv("../../data/sample_info/sample_level_covariates_atac_with_frag.tsv", sep = "\t")
sample_info <- sample_info %>% filter(coarse_cluster_name == 'fusion.Type_1')
# current list of samples - 284
sample_list <- sample_info$SNG.1ST
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
## get L1 and L2
L_1= dat_281 %>% select(hit1)
L_2= dat_281 %>% select(hit2)

# --------------------peak info-----------------------
nfr_peak_id = nfr_peak
nuc_peak_id = nuc_peak
## count mat locations
nfr_mat=paste0(opts$count_mat_dir, celltype, '_nfr_normalized.rds')
nuc_mat=paste0(opts$count_mat_dir, celltype, '_nuc_normalized.rds')
#HMMRATAC_nfr_mat = "../../data/normalized_counts/HMMRATAC_NFR_counts_multiotsu_normalized.rds"
#HMMRATAC_nuc_mat = "../../data/normalized_counts/HMMRATAC_Nucleosomal_counts_multiotsu_normalized.rds"
#NucleoATAC_nfr_mat = "../../data/normalized_counts/NuclaoATAC_NFR_counts_multiotsu_byreads_normalized.rds"
#NucleoATAC_nuc_mat = "../../data/normalized_counts/NucleoATAC_Nucleosomal_counts_multiotsu_200ext_04282023_filtered_normalized.rds"
# if (trait == "HMMRATAC"){
#         nfr_count_mat = readRDS(HMMRATAC_nfr_mat)
#         nuc_count_mat = readRDS(HMMRATAC_nuc_mat)
# } else if (trait == "NucleoATAC"){
#         nfr_count_mat = readRDS(NucleoATAC_nfr_mat)
#         nuc_count_mat = readRDS(NucleoATAC_nuc_mat)
# }
# read in
nfr_count_mat = readRDS(nfr_mat)
nuc_count_mat = readRDS(nuc_mat)
# get data points
nfr_col = nfr_count_mat %>% as.data.frame() %>% select(all_of(nfr_peak_id))
nuc_col = nuc_count_mat %>% as.data.frame() %>% select(all_of(nuc_peak_id))
# get sample id
nfr_names = nfr_col %>% rownames()
nfr_names = sapply(strsplit(nfr_names, "__"), function(x) x[[1]])
nuc_names = nuc_col %>% rownames()
nuc_names = sapply(strsplit(nuc_names, "__"), function(x) x[[1]])
# add names col
rownames(nfr_col) = nfr_names
rownames(nuc_col) = nuc_names
peak_df = merge(nfr_col, nuc_col, by = 'row.names', all = TRUE)
# rename cols
colnames(peak_df) = c("Sample", "NFR_peak", "Nuc_peak")

# --------------------Covariates info-----------------------
cov_file = paste0(cov_file, celltype, '_nfr_FPC1_GPC5.txt') # either nfr or nuc should be fine
cov_df = read.table(cov_file, , header=TRUE, check.names=F)
#load(cov_file) # should be C
C = cov_df %>% filter(id != 'FPC0')
# remove the names of covs
C_mat = C[,-1]
# transpose it to be nxp
C_mat = t(C_mat)
#C_mat %>% head()
# order the samples in the matrix
C_mat_ordered = merge(L_1, C_mat, by.x=0, by.y=0) %>% select(-c(1,2))
# the same C mat should be applied to all hypothesis
C = C_mat_ordered

# --------------------Assemble cit df-----------------------
## one driven by nfr_hit - L1
cit_df_1 = merge(L_1, peak_df, by.x = 0, by.y = "Sample")
## one driven by nuc_hit - L2
cit_df_2 = merge(L_2, peak_df, by.x = 0, by.y = "Sample")

# --------------------Assemble 4 hypothesis-----------------------
# H1: hit1 -> nfr_peak -> nuc_peak
L1 = cit_df_1 %>% select(hit1) %>% as.matrix()
G1 = cit_df_1 %>% select("NFR_peak") %>% as.matrix()
T1 = cit_df_1 %>% select("Nuc_peak") %>% as.matrix()
h1_results=cit.cp(L1,G1,T1,C, n.perm=100)

# H2: hit1 -> nuc_peak -> nfr_peak
G2 = cit_df_1 %>% select("Nuc_peak") %>% as.matrix()
T2 = cit_df_1 %>% select("NFR_peak") %>% as.matrix()
h2_results=cit.cp(L1,G2,T2,C, n.perm=100)

# H3: hit2 -> nfr_peak -> nuc_peak
L2 = cit_df_2 %>% select(hit2) %>% as.matrix()
G3 = cit_df_2 %>% select("NFR_peak") %>% as.matrix()
T3 = cit_df_2 %>% select("Nuc_peak") %>% as.matrix()
h3_results=cit.cp(L2,G3,T3,C, n.perm=100)

# H4: hit2 -> nuc_peak -> nfr_peak
G4 = cit_df_2 %>% select("Nuc_peak") %>% as.matrix()
T4 = cit_df_2 %>% select("NFR_peak") %>% as.matrix()
h4_results=cit.cp(L2,G4,T4,C, n.perm=100)

# --------------------Save results-----------------------
results = list()
results$H1 = h1_results
results$H2 = h2_results
results$H3 = h3_results
results$H4 = h4_results
# save dat list
#ave(results, file = paste0(outdir, 'PPH4_threshold', PPH4_threshold, '/', cit_pair, '_cit_perm_nocov.Rda'))
save(results, file = paste0(outdir, 'cit_PPH4_threshold', PPH4_threshold, '/', cit_pair, '_cit_perm_cov.Rda'))
