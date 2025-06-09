#!/usr/bin/env Rscript

# libraries
library(limma)
library(edgeR)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggfortify)
library(optparse)
library(RNOmni)
options(stringsAsFactors=FALSE)

set.seed(1) 

option_list <- list(
  make_option(c("--counts"), type = "character", help = "featureCounts counts file."),
  make_option(c("--counts_norm"), type = "character", help = "Output rds file."),
  make_option(c("--outdirQC"), type = "character", help = "Location to store output QC plots."),
  make_option(c("--outdirBED"), type = "character", help = "Location to store output peak bed files."),
  make_option(c("--sample_info_file"), type = "character", help = "Sample information."),
  make_option(c("--celltype"), type = "character", help = "Which celltype"),
  make_option(c("--region"), type = "character", help = "which region? nfr/nuc")
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]",
                              option_list = option_list, add_help_option = TRUE)
opts <- parse_args(option_parser)

# read in parameters
# test opts
#opts$counts = "../../results/results_02062024/counts_file_filtered/Type_2a_nfr_counts_filtered.txt"
#opts$sample_info_file = "../../data/sample_info/sample_level_covariates_atac.tsv"
#opts$outdir = "../../results/results_test/"
#opts$counts_norm = "../../results/results_test/Type_2a_nfr_counts_nomalized.rds"

# read data
count_matrix_file = opts$counts
sample_info_file = opts$sample_info_file

cat("Reading count matrix:", count_matrix_file, "\n")
countMat = read.csv(count_matrix_file, sep = '\t', check.names=F)

cat("Reading sample info file:", sample_info_file, "\n")
sampleDF = read.table(sample_info_file, sep='\t', header = TRUE)
# filter and match for the samplDF
sampleDF <- sampleDF[c('modality', 'batch', 'SAMPLE')]
sampleDF <- sampleDF[!duplicated(sampleDF),]
# order sample
order <- vapply(strsplit(colnames(countMat), "__"), '[', 1, FUN.VALUE = character(1))
order <- order[-1]
# remove unqualified samples - deindentified samples
order <- order[order != "donor1"]
order <- order[order != "donor2"]
order <- order[order != "donor3"]
sampleDF <- sampleDF[match(order, sampleDF$SAMPLE),]
cat("Extracting peak info\n")
formatPeakName = function(x) {
  y = unlist(strsplit(x, split=".", fixed=TRUE))
  n = length(y)
  start = as.numeric(y[(n-1)])
  end = as.numeric(y[n])
  chr = paste(y[1:(n-2)],collapse="_")
  peakname = paste0(chr, ":", as.character(start), "-", as.character(end))
  return(data.frame(peakname, chr, start, end))
}
# tidy the countmat
## rownames
countMat_tidy <- countMat[,-1]
rownames(countMat_tidy) <- countMat[,1]
## colnames
correct_names <- lapply(colnames(countMat_tidy), function(x) strsplit(x, split='__')[[1]][1])
colnames(countMat_tidy) <- correct_names
# run the peakformating function
peakList =  lapply(rownames(countMat_tidy), formatPeakName)
peakDF = do.call(rbind, peakList)
peakDF$length = peakDF$end - peakDF$start
row.names(peakDF) = rownames(countMat_tidy)

cat("Filtering low peaks using TPM\n")
### Calculate TPM (note: this is from https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/)
countToTpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

tpm = apply(countMat_tidy, MARGIN=2, function(x) countToTpm(x, effLen=peakDF[row.names(countMat_tidy),"length"]))
keep_peaks = row.names(tpm)[apply(tpm>1, MARGIN=1, FUN=mean) > 0.1]
countMat_filtered = countMat_tidy[keep_peaks,]

cat("Removing unqualified samples\n")
countMat_filtered <- countMat_filtered %>% select(-c('12004', '22011', '32071'))

cat("Normalizing count matrix using TMM\n")

### Create edgeR object
y_raw <- DGEList(counts=countMat_filtered, samples=sampleDF)
lcpm_raw <- cpm(y_raw, log=TRUE)

### Calculate log-CPM values before and after TMM normalization
y_norm <- calcNormFactors(y_raw, method="TMM")
cpm_norm <- cpm(y_norm, log=FALSE)
lcpm_norm <- cpm(y_norm, log=TRUE)

cat("Inverse rank normalizing count matrix \n")
inverseNormalize = function(x) {
  y = qnorm((rank(x, na.last="keep")-0.5) / sum(!is.na(x)))
  return(y)
}
#cpm_inversenorm = apply(cpm_norm, MARGIN=1, FUN=inverseNormalize)
# randomly break the tie
cpm_inversenorm = apply(cpm_norm, MARGIN=1, function(x) FUN=RankNorm(x,ties.method='random'))

cat("TMM-normalize and inverse normalized autosomal count matrix: ", opts$counts_norm, "\n")
saveRDS(cpm_inversenorm, file=opts$counts_norm)

cat("Plot QCs \n")
pca_df <- merge(cpm_inversenorm, sampleDF, by.x='row.names', by.y='SAMPLE')
pca_res <- pca_df %>% select(-c('Row.names', 'modality', 'batch')) %>% prcomp(scale. = TRUE)
# PC1vs2
g1 <- autoplot(pca_res, data = pca_df, colour = 'batch')
# PC2vs3
g2 <- autoplot(pca_res, data = pca_df, colour = 'batch', x=2, y=3)
# PC3vs4
g3 <- autoplot(pca_res, data = pca_df, colour = 'batch', x=3, y=4)
# variance plot
png(file=paste0(opts$outdirQC, 'pca_cumvar.png'))
summary(pca_res)$importance[3,] %>% plot(xlab='PCs', ylab='Cumulative Proportion of Variance Explained')
dev.off()

ggsave(filename=paste0(opts$outdirQC, 'pc1vs2.png'), g1, device='png')
ggsave(filename=paste0(opts$outdirQC, 'pc2vs3.png'), g2, device='png')
ggsave(filename=paste0(opts$outdirQC, 'pc3vs4.png'), g3, device='png')

# save peak df
cat(paste0('track name=', opts$celltype, '_', opts$region, 'Peaks description=\"', opts$celltype, ' ', 
            opts$region, 'Peaks NucleoATAC"'),
            file=paste0(opts$outdirBED, opts$celltype, '_', opts$region, '.bed'))

peak_bed <- peakDF %>% select(c('chr', 'start', 'end'))
cat(" \n",file=paste0(opts$outdirBED, opts$celltype, '_', opts$region, '.bed'), append=TRUE)
options(scipen=999)
write.table(peak_bed, file=paste0(opts$outdirBED, opts$celltype, '_', opts$region, '.bed'),
            sep="\t",  col.names=FALSE, row.names=FALSE, quote=FALSE, append=TRUE)