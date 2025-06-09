#!/usr/bin/env Rscript

# -----------libraries-----------
options(stringsAsFactors=FALSE)
library(susieR)
library(dplyr)
library(vcfR)
library(coloc)
library(Rfast)
library(optparse)

option_list <- list(
    make_option(c("--GWAS_annot"), type = "character", help = "[Required] Which GWAS summary stats?"),
    make_option(c("--gwas_file_path"), type = "character", help = "[Required] Building susie models for which GWAS signals"),
    make_option(c("--LD_file"), type = "character", help = "[Required] Lead variant output filename"),
    make_option(c("--ukbb_path"), type = "character", help = "[Required] Using which VCF information, ukbb?"),
    make_option(c("--min_corr"), type = "numeric", help = "[Required] min abs corr between any snp pairs"),
    make_option(c("--num_L"), type = "numeric", help = "[Required] Max number of LD vars"),
    make_option(c("--max_it"), type = "numeric", help = "[Required] Max number of susie iterations"),
    make_option(c("--outdir"), type = "character", help = "[Required] Susie models output directory")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# variables
GWAS_annot = opts$GWAS_annot
LD_file = opts$LD_file
ukbb_path = opts$ukbb_path
gwas_file_path = opts$gwas_file_path
num_L = opts$num_L
min_corr = opts$min_corr
max_it= opts$max_it
outdir = opts$outdir

LD_file_temp = strsplit(LD_file, split="gwas_")[[1]][2]
chr = strsplit(LD_file_temp, split="-")[[1]][1]

# read variables in
#gwas_leadvar = read.table(gwas_leadvar_path, sep = "\t", header=TRUE)

# filter to chr
#gwas_leadvar_chr = gwas_leadvar %>% filter(Chr == strsplit(chr, split='chr')[[1]][2])


# ukbb vcf file
ukbb_vcf_path = paste0(ukbb_path, chr, ".imputed.poly.vcf.gz")

# gwas summary stats
gwas_stats = data.table::fread(gwas_file_path, header=TRUE)

# filter gwas summary stats to chr
gwas_stats_chr = gwas_stats %>% filter(`#snp_chrom` == chr)



print("-------------Assembling genotype information-------------")
# creating sub vcfs
LD_chr = chr
LD_start = strsplit(LD_file_temp, split="-")[[1]][2]
LD_end = strsplit(LD_file_temp, split="-")[[1]][3]
LD_pos = paste0(LD_chr, ":", LD_start, "-", LD_end)
#LD_file = paste0("gwas_", LD_chr, "-", LD_start, "-", LD_end)
command = paste0("tabix -h ", ukbb_vcf_path, " ")
command = paste0(command, LD_pos)
# change the directory to "../../data" in smk!!!
command = paste0(command, " | bgzip > ../../data/sample_info/temp_vcfs/", LD_file, ".vcf.gz")
print(command)
system(command, intern=TRUE)
print("Done with running tabix")
# get the genotype inforamtion from sub vcfs
sub_vcf_path = paste0("../../data/sample_info/temp_vcfs/", LD_file, ".vcf.gz")
# read vcf
vcf = read.vcfR(sub_vcf_path)
print("Done with reading vcf")
# get id pos and allele info
# from the fix region
vcf_info_df = getFIX(vcf) %>% as.data.frame()
# col 1
vcf_info_df$pos_id_ref_alt = apply(vcf_info_df, 1, FUN=function(x) paste0(x['POS'],"_",x['ID'],"_",x['REF'],"_",x['ALT']))
# col 2
vcf_info_df$pos_id_alt_ref = apply(vcf_info_df, 1, FUN=function(x) paste0(x['POS'],"_",x['ID'],"_",x['ALT'],"_",x['REF']))
dosages = extract.gt(vcf, element="DS", IDtoRowNames=FALSE)
#dat <- mutate_all(as.data.frame(dosages, row.names=vcf_info_df$pos_id_ref_alt), function(x) as.numeric(x))
rownames_nodup = ave(as.character(vcf_info_df$pos_id_ref_alt), vcf_info_df$pos_id_ref_alt, FUN=function(x) if (length(x)>1) paste0(x[1], '(', seq_along(x), ')') else x[1])
dat <- mutate_all(as.data.frame(dosages, row.names=rownames_nodup), function(x) as.numeric(x))
# transpose dataframe
dat_gwas <- t(dat)
# turn into a dataframe
dat_gwas <- as.data.frame(dat_gwas)

print("---------Assembling coloc datalist: GWAS part---------")

# snps in the dosage df
snp_list = vcf %>% getID()
print(paste("Number of snps in the UKBB vcf:", length(snp_list)))
coloc_gwas_df = gwas_stats_chr %>% filter(SNP %in% snp_list)
print(paste("Number of snps overlap between UKBB vcf and GWAS stats:", nrow(coloc_gwas_df)))
# Fix the alignment
print("---------Fixing GWAS vcf alignment---------")
# col 1
coloc_gwas_df$pos_id_ref_alt = apply(coloc_gwas_df, 1, FUN=function(x) gsub(" ", "", paste0(x['snp_end'],"_",x['SNP'],"_",x['NEA'],"_",x['EA'])))
# col 2
coloc_gwas_df$pos_id_alt_ref = apply(coloc_gwas_df, 1, FUN=function(x) gsub(" ", "", paste0(x['snp_end'],"_",x['SNP'],"_",x['EA'],"_",x['NEA'])))


# init track lists
beta_modified = c()
keep_list = c()

# looping through the coloc_gwas_df is better option
# to check the ref/alt vs. NEA/EA alignment
for (j in 1: nrow(coloc_gwas_df)){
    # rsid col
    rsid = coloc_gwas_df[j,]$SNP
    # full id col
    fullid = coloc_gwas_df[j,]$pos_id_ref_alt
    # remove the space in fullid
    #fullid = gsub(" ", "", fullid)
    # match the vcf info
    vcf_snp_df = vcf_info_df %>%
                filter(ID == rsid)
    # keep the designated rows
    if (nrow(vcf_snp_df) > 1){
        print(rsid)
        #print(fullid)
        #print(nrow(vcf_snp_df))
        vcf_snp_df = vcf_snp_df %>% filter(pos_id_ref_alt == fullid | pos_id_alt_ref == fullid)
        #print(nrow(vcf_snp_df))
    }
    
    # if this row is a match
    if (fullid == vcf_snp_df$pos_id_ref_alt){
        # keep the beta
        beta_modified = c(beta_modified, coloc_gwas_df[j,]$beta)
        # keep the snp
        keep_list = c(keep_list, 1)
    }
    else if (fullid == vcf_snp_df$pos_id_alt_ref){
        # keep the beta
        beta_modified = c(beta_modified, coloc_gwas_df[j,]$beta*(-1.0))
        # keep the snp
        keep_list = c(keep_list, 1)
    } else{
        beta_modified = c(beta_modified, coloc_gwas_df[j,]$beta)
        keep_list = c(keep_list, 0)
    }
}
coloc_gwas_df$beta_modified = beta_modified
coloc_gwas_df$keep = keep_list
print(paste("Number of snps in GWAS stats before checking alignment:", nrow(coloc_gwas_df)))
# drop those without alignment
coloc_gwas_df = coloc_gwas_df %>% filter(keep == 1)
print(paste("Number of snps in GWAS stats after checking alignment:", nrow(coloc_gwas_df)))
# Filter the vcf dosage information
print("---------Filtering for snps in dosage df---------")
# union of potential ids to keep
keep_fullids = union(coloc_gwas_df$pos_id_ref_alt, coloc_gwas_df$pos_id_alt_ref)
LD_dat = dat_gwas[, which((names(dat_gwas) %in% keep_fullids)==TRUE)]
print("---------Sanity check---------")
if (length(names(LD_dat)) != nrow(coloc_gwas_df)){
    print("Warning: potential misalignment in the VCF dosgae inforamtion and GWAS summary stats")
    print("Please check the pos, rsid, ref, and alt alleles' inforation")
} else{print("PASS!")}
# Rename the dosage df
dosage_snps = read.table(text=names(LD_dat),sep='_',strip.white = TRUE, col.names=c("Pos", "SNP","ref", "alt"))$SNP
colnames(LD_dat) = dosage_snps
print("-------------Generating the LD mat -------------")
# calculate cor
LD_dat_m <- cor(LD_dat)
print("-------------Assembling colco datalist -------------")
coloc_d1 <- list()
coloc_d1$beta = coloc_gwas_df$beta_modified
coloc_d1$varbeta = (coloc_gwas_df$se)^2
coloc_d1$snp = coloc_gwas_df$SNP
coloc_d1$position = coloc_gwas_df$snp_end
coloc_d1$type = "cc"
coloc_d1$LD = LD_dat_m
coloc_d1$pvalues = coloc_gwas_df$pval
# save coloc data
save(coloc_d1, file=paste0(outdir, "/", GWAS_annot, "_susie_models/", LD_file, "_susie_dat.Rda"))
# Run susie
S1 <- susie_rss(
    bhat = coloc_d1$beta,
    shat = coloc_d1$varbeta,
    R = coloc_d1$LD,
    estimate_residual_variance = FALSE,
    L = num_L,
    min_abs_corr = min_corr,
    max_iter=max_it)
# Save susie model
save(S1, file=paste0(outdir, "/", GWAS_annot, "_susie_models/", LD_file, ".susie.Rda"))
pipv = susie_get_pip(S1)
pip_df = data.frame(list("pip"=pipv, "snp"=coloc_d1$snp))
write.table(pip_df, file=paste0(outdir, "/", GWAS_annot, "_susie_models/", LD_file, ".susie_pipdf.tsv"),
            sep='\t', quote=F)
# is there 95% credible set?
## how many 95% csets?
n95 = length(S1$sets$cs)
## best set
best = row.names(S1$sets$purity[which.max(S1$sets$purity$min.abs.corr),])
## size of the main set
len_1_95 = length(S1$sets$cs[[best]])
## save the main 95% cset
write.table(pip_df[S1$sets$cs[[best]],], 
            file=paste0(outdir, "/", GWAS_annot, "_susie_models/", LD_file, '.susie_cset95.tsv'), 
            sep='\t', quote=F)
print(paste(n95, "95% credible sets found. Length of the first 95% credible set =", len_1_95))
print("-------------DONE for this lead var -------------")
