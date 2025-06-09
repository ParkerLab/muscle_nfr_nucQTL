library(dplyr)
library(tidyr)
library(stringr)
library(seqinr)
library(vcfR)
library(optparse)

option_list <- list(
    make_option(c("--inpeak_dir"), type = "character", help = "[Required] Inpeak snps results dir"),
    make_option(c("--var_info"), type = "character", help = "[Required] Varaints' ref/alt information"),
    make_option(c("--hg38_fasta"), type = "character", help = "[Required] hg38 fasta file"),
    make_option(c("--part1_bed"), type = "character", help = "[Required] bed file before snp 75bp"),
    make_option(c("--part2_bed"), type = "character", help = "[Required] bed file after snp 75bp"),
    make_option(c("--part1_fasta"), type = "character", help = "[Required] fasta file before snp 75bp"),
    make_option(c("--part2_fasta"), type = "character", help = "[Required] fasta file after snp 75bp"),
    make_option(c("--output_open"), type = "character", help = "[Required] open nfrQTLs FASTA output dir"),
    make_option(c("--output_close"), type = "character", help = "[Required] close nfrQTLs FASTA output dir")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# input
inpeak_file = opts$inpeak_dir
var_file = opts$var_info
hg38 = opts$hg38_fasta
# output
part1_bed = opts$part1_bed
part2_bed = opts$part2_bed
part1_fasta = opts$part1_fasta
part2_fasta = opts$part2_fasta
openFasta <- opts$output_open
closeFasta <- opts$output_close

# functions
# read the fasta files
read_fasta <- function(raw_file){
    # read in
    data_raw = read.table(raw_file)
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
    #combined_data <- combined_data %>% unique()
    return(combined_data)
}

# read in files
# a list of inpeak qtls
df = read.table(inpeak_file) %>%
                filter(in_peak == 1) %>% unique()
# Filter rows with the same 'phe_id'
dup_phe_ids <- df$phe_id[duplicated(df$phe_id) | duplicated(df$phe_id, fromLast = TRUE)]
dup_df <- df[df$phe_id %in% dup_phe_ids, ]

# Find the row with the highest 'pip' for each 'phe_id', choosing the first one in case of ties
highest_pip_df <- dup_df[order(-dup_df$pip), ]
highest_pip_df <- highest_pip_df[!duplicated(highest_pip_df$phe_id), ]

nodup_df <- df[!df$phe_id %in% dup_phe_ids, ]

inpeak_df = rbind(highest_pip_df, nodup_df)
# remove those that are indels
inpeak_df = inpeak_df %>% filter(var_from == var_to)

# varinfo
variant_data <- read.table(var_file)
# remove indels
variant_data <- variant_data %>% filter((str_length(Ref) == 1) & (str_length(Alt) == 1))
#variant_data

# get the columns
inpeak_var_pos <- inpeak_df %>% select(c('var_chr', 'var_from', 'var_to'))

# get the partone df -> bed
part1_df <- inpeak_var_pos
part1_df$chrom <- inpeak_var_pos$var_chr
part1_df$chromStart <- inpeak_var_pos$var_from - 76
part1_df$chromEnd <-inpeak_var_pos$var_from -1
# save to bedfile
part1_df <- part1_df %>% select(c('chrom', 'chromStart', 'chromEnd'))
write.table(part1_df, part1_bed, sep = '\t', quote = F, row.names = F, col.names = F)

# get the parttwo df -> bed
part2_df <- inpeak_var_pos
part2_df$chrom <- inpeak_var_pos$var_chr
part2_df$chromStart <- inpeak_var_pos$var_from
part2_df$chromEnd <-inpeak_var_pos$var_from + 75
# save to bedfile
part2_df <- part2_df %>% select(c('chrom', 'chromStart', 'chromEnd'))
write.table(part2_df, part2_bed, sep = '\t', quote = F, row.names = F, col.names = F)
# get fasta
# Construct the bedtools command1
command1 <- paste("bedtools getfasta",
                 "-fi", hg38,
                 "-bed", part1_bed,
                 "-fo", part1_fasta)
system(command1)
# get fasta
# Construct the bedtools command2
command2 <- paste("bedtools getfasta",
                 "-fi", hg38,
                 "-bed", part2_bed,
                 "-fo", part2_fasta)
system(command2)
# read in the seq
fasta1_df <- read_fasta(part1_fasta)
fasta2_df <- read_fasta(part2_fasta)
# check if they are of the same length
if (!nrow(fasta1_df) == nrow(fasta2_df)) {
  stop("Error: FASTA files have different length")
}
if (!nrow(fasta1_df) == nrow(inpeak_df)) {
  stop("Error: FASTA file and inpeak information have different length")
}

# Initialize (overwrite) files
cat("", file=openFasta)
cat("", file=closeFasta)
for(i in 1:nrow(inpeak_df)){
#for(i in 1:1){
    print(i)
    # get the entry row
    var_entry = variant_data %>% filter(Chr == inpeak_df[i,]$var_chr,
                        Pos == inpeak_df[i,]$var_from)
    # get the vcf ref
    vcf_ref = var_entry$Ref
    # get the vcf alt
    vcf_alt = var_entry$Alt
    # assemble seqs
    mod_seq = paste0(fasta1_df[i,]$Seq, vcf_alt, fasta2_df[i,]$Seq)
    ori_seq = paste0(fasta1_df[i,]$Seq, vcf_ref, fasta2_df[i,]$Seq)
    # generate new peak name
    peakname = paste0('>', part1_df[i,]$chrom, ':', part1_df[i,]$chromStart, '-', part2_df[i,]$chromEnd)
    # decide if this is open allele or close
    if (inpeak_df[i,]$slope > 0){
        print('open')
        # open
        write(peakname, file=openFasta, append=TRUE, sep='\n')
        write(mod_seq, file=openFasta, append=TRUE, sep='\n')
        write(peakname, file=closeFasta, append=TRUE, sep='\n')
        write(ori_seq, file=closeFasta, append=TRUE, sep='\n')
    } else {
        print('close')
        # close
        write(peakname, file=openFasta, append=TRUE, sep='\n')
        write(ori_seq, file=openFasta, append=TRUE, sep='\n')
        write(peakname, file=closeFasta, append=TRUE, sep='\n')
        write(mod_seq, file=closeFasta, append=TRUE, sep='\n')
    }
}    