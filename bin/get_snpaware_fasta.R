library(dplyr)
library(tidyr)
library(stringr)
library(seqinr)
library(vcfR)

library(optparse)

option_list <- list(
    make_option(c("--inpeak_dir"), type = "character", help = "[Required] Inpeak snps results dir"),
    make_option(c("--var_info"), type = "character", help = "[Required] Varaints' ref/alt information"),
    make_option(c("--nfrqtl_fasta"), type = "character", help = "[Required] coloc nfrQTLs' fasta file"),
    make_option(c("--output_open"), type = "character", help = "[Required] open nfrQTLs FASTA output dir"),
    make_option(c("--output_close"), type = "character", help = "[Required] close nfrQTLs FASTA output dir")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

openFasta <- opts$output_open
closeFasta <- opts$output_close

# read in files
# a list of inpeak qtls
df = read.table(opts$inpeak_dir) %>%
                filter(in_peak == 1) %>% unique()
# variant information
variant_data <- read.table(opts$var_info)
# remove indels
variant_data <- variant_data %>% filter((str_length(Ref) == 1) & (str_length(Alt) == 1))
#variant_data
# Filter rows with the same 'phe_id'
dup_phe_ids <- df$phe_id[duplicated(df$phe_id) | duplicated(df$phe_id, fromLast = TRUE)]
dup_df <- df[df$phe_id %in% dup_phe_ids, ]

# Find the row with the highest 'pip' for each 'phe_id', choosing the first one in case of ties
highest_pip_df <- dup_df[order(-dup_df$pip), ]
highest_pip_df <- highest_pip_df[!duplicated(highest_pip_df$phe_id), ]

nodup_df <- df[!df$phe_id %in% dup_phe_ids, ]
# merge two parts
inpeak_df = rbind(highest_pip_df, nodup_df)

# raw per chr qtl fasta file
data_raw=read.table(opts$nfrqtl_fasta)
# Separate the odd and even rows
odd_rows <- data_raw[seq(1, nrow(data_raw), by = 2), ]
even_rows <- data_raw[seq(2, nrow(data_raw), by = 2), ]

# Combine the odd and even rows into a new data frame
combined_data <- data.frame(Peak = odd_rows, Seq = even_rows)

# match the qtl phe_id
combined_data$Peak_name = gsub(">", "", combined_data$Peak)
# reformatting
combined_data$Peak_name = gsub(":", ".", combined_data$Peak_name)
combined_data$Peak_name = gsub("-", ".", combined_data$Peak_name)
# drop duplicates
combined_data <- combined_data %>% unique()

# Initialize (overwrite) files
cat("", file=openFasta)
cat("", file=closeFasta)
#for(i in 1:nrow(inpeak_df)){
for(i in 1:1:nrow(inpeak_df)){
    # get the entry row
    var_entry = variant_data %>% filter(Chr == inpeak_df[i,]$var_chr,
                        Pos == inpeak_df[i,]$var_from)
    # get the entry from fasta side
    entry = combined_data %>% 
            filter(Peak_name == paste0(inpeak_df[i,]$phe_chr, '.',
                                       inpeak_df[i,]$phe_from, '.',
                                       inpeak_df[i,]$phe_to))
    if (nrow(var_entry) == 0 | nrow(entry) == 0){
        next
    }
    # check the if ref allele match
    fasta_ref = str_split(entry$Seq, '')[[1]][(inpeak_df[i,]$var_from-inpeak_df[i,]$phe_from):(inpeak_df[i,]$var_from-inpeak_df[i,]$phe_from)]
    # get the vcf ref
    vcf_ref = var_entry$Ref
    # get the vcf alt
    vcf_alt = var_entry$Alt

    if (inpeak_df[i,]$var_from-inpeak_df[i,]$phe_from == 0){
        mod_seq = paste0(vcf_alt, entry$Seq)
        ori_seq = paste0(vcf_ref, entry$Seq)
    } else {
        # modify the Seq
        ## get the old part 1
        part1 = str_split(entry$Seq, '')[[1]][0:(inpeak_df[i,]$var_from-inpeak_df[i,]$phe_from-1)]
        part1 = paste(part1, collapse='')
        ## get the old seq part 2
        part2 = str_split(entry$Seq, '')[[1]][(inpeak_df[i,]$var_from-inpeak_df[i,]$phe_from+1):151]
        part2 = paste(part2, collapse='')
        ## get modified seq
        mod_seq = paste0(part1, vcf_alt, part2)
        # remove NA just in case
        mod_seq = gsub('NA', '', mod_seq)
        ori_seq = entry$Seq
    }

    # decide if this is open allele or close
    if (inpeak_df[i,]$slope > 0){
        print('open')
        # open
        write(entry$Peak, file=openFasta, append=TRUE, sep='\n')
        write(mod_seq, file=openFasta, append=TRUE, sep='\n')
        write(entry$Peak, file=closeFasta, append=TRUE, sep='\n')
        write(ori_seq, file=closeFasta, append=TRUE, sep='\n')
    } else {
        print('close')
        # close
        write(entry$Peak, file=openFasta, append=TRUE, sep='\n')
        write(ori_seq, file=openFasta, append=TRUE, sep='\n')
        write(entry$Peak, file=closeFasta, append=TRUE, sep='\n')
        write(mod_seq, file=closeFasta, append=TRUE, sep='\n')
    }
    
    
}