#!/usr/bin/env Rscript

# Description: This script uses motifbreakR to analyze motifs in genomic data.
# Author: Alice Wang
# Date: 05012025
# Usage: Rscript motif_mb.R [arguments]

library(VariantAnnotation)
library(GenomicRanges)
library(Biostrings)
library(dplyr)
library(vcfR)
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
library(MotifDb)
library(BSgenome)
library(universalmotif)
library(TFBSTools)
library(S4Vectors)
library(stringr)
library(optparse)

option_list <- list(
    make_option(c("--inpeak_df"), type = "character", help = "[Required] Inpeak snps results dir"),
    make_option(c("--output"), type = "character", help = "[Required] output motifbreakR results dir"),
    make_option(c("--motif_db"), type = "character", help = "[Required] motif database")
)
option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T)
opts = parse_args(option_parser)

# Path to your log-odds PWM file
file_path <- opts$motif_db
# Read all lines
lines <- readLines(file_path)
motif_lines <- which(grepl("^>", lines))

# Prepare containers
motif_list <- list()
providerName <- character()
providerId <- character()
geneSymbol <- character()
geneId <- character()
proteinId <- character()
bindingSequence <- character()
subtype <- character()
source <- character()
quality <- character()

# Iterate through each motif block
for (i in seq_along(motif_lines)) {
  start <- motif_lines[i]
  end <- if (i == length(motif_lines)) length(lines) else motif_lines[i + 1] - 1

  # Extract motif name
  motif_name <- sub("^>", "", lines[start])
  id <- paste0("custom-", motif_name)

  # Read matrix rows (assumes order A, C, G, T)
  pwm_matrix <- do.call(rbind, lapply(lines[(start + 1):end], function(line) {
    as.numeric(strsplit(line, "\\s+")[[1]])
  }))
  colnames(pwm_matrix) <- c("A", "C", "G", "T")
  
  # Store matrix
  motif_list[[id]] <- t(pwm_matrix)
  

  # Fill dummy metadata
  providerName <- c(providerName, str_split(motif_name, '\\.', 5)[[1]][2])
  providerId   <- c(providerId, motif_name)
  geneSymbol   <- c(geneSymbol, str_split(motif_name, '\\.', 5)[[1]][1])
  geneId       <- c(geneId, NA_character_)
  proteinId    <- c(proteinId, NA_character_)
  bindingSequence <- c(bindingSequence, paste(rep("N", ncol(pwm_matrix)), collapse = ""))
  subtype <- c(subtype, str_split(motif_name, '\\.', 5)[[1]][3])
  source <- c(source, str_split(motif_name, '\\.', 5)[[1]][4])
  quality <- c(quality, str_split(motif_name, '\\.', 5)[[1]][5])
}

# Construct metadata
listData <- DataFrame(
  providerName = providerName,
  providerId = providerId,
  dataSource = source,
  geneSymbol = geneSymbol,
  geneId = geneId,
  geneIdType = NA_character_,
  proteinId = proteinId,
  proteinIdType = NA_character_,
  organism = "Hsapiens",
  sequenceCount = NA_integer_,
  bindingSequence = bindingSequence,
  bindingDomain = NA,
  tfFamily = NA,
  experimentType = subtype,
  pubmedID = NA_character_,
  subtype = subtype,
  quality = quality
)
rownames(listData) <- names(motif_list)

# Create MotifDb object using internal constructor
MotifList <- getFromNamespace("MotifList", "MotifDb")
motifdb_obj <- MotifList(motif_list, listData)

# start analysis
# read variatnts
var_df = read.table(opts$inpeak_df, header = TRUE, sep = "\t")
var_df$snp[var_df$snp == 'chr1:21937420:G:C'] <- 'rs192259893'

variants <- snps.from.rsid(rsid = var_df$snp,
               dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38,
               search.genome = BSgenome.Hsapiens.UCSC.hg38)


motifbreakr.results <- motifbreakR(snpList = variants, 
                                   filterp = TRUE,
                                   pwmList = motifdb_obj,
                                   threshold = 0.1,
                                   method = "default",
                                   bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                   BPPARAM = BiocParallel::SerialParam())

print("Done")
# filter results
results = motifbreakr.results[motifbreakr.results$effect == 'strong']
results = results[grepl(".0.", results$providerId, ignore.case = TRUE)]
#results = results[results$SNP_id == 'rs646527']
#results = calculatePvalue(results, BPPARAM = BiocParallel::SerialParam())

saveRDS(results, file = opts$output)