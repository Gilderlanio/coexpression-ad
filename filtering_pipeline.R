#!/usr/bin/env Rscript
#'
#' This is the R Script to perform filtering step for co-expression analysis.
#'
#' Corresponding authors: Prof. Dr. Gilderlanio Santana de Araujo (gilderlanio [at] gmail.com)
#'                        Ms.C Arthur Ribeiro dos Santos (arthurrdsantos [at] outlook.com)

# Load libraries
library(stringr)
library(hgnc)

# Load data
AD_GSE125583 <- read.delim2("~/counts/raw_data/AD_GSE125583_raw.txt",
                            header = T, row.names = 1,
                            sep = "\t", stringsAsFactors = F)

NC_GSE125583 <- read.delim2("~/counts/raw_data/CN_GSE125583_raw.txt",
                            header = T, row.names = 1,
                            sep = "\t", stringsAsFactors = F)

# Create NC+AD data
All_GSE125583 <- cbind(AD_GSE125583, NC_GSE125583)


# Select group for the process

## this will be repeated later for NC_GSE125583 and All_GSE125583
input_data <- AD_GSE125583


# Filter by ensemblID match

## get ensemblID names
ensemblID <- row.names(input_data)
ensemblID <- t(as.data.frame(str_split(ensemblID, pattern = "[.]")))
ensemblID <- ensemblID[,1]

## check for hgnc_complete_set.txt file and download if necessary
if(!file.exists("~/hgnc_data/hgnc_complete_set.txt")){
  download_archive()
}

## load HGNC file
hgnc_full_list = read.delim2("~/hgnc_data/hgnc_complete_set.txt",
                             header = T, sep = "\t", stringsAsFactors = F)

## check for matches in HGNC file
k <- ensemblID %in% hgnc_full_list$ensembl_gene_id

## generate dataset filtered by ensemblID match
GSE125583_Symbol <- GSE125583[k,]

## renaming rows to get gene Symbol
k <- na.omit(match(ensemblID, hgnc_full_list$ensembl_gene_id))
row.names(GSE125583_Symbol) <- hgnc_full_list$symbol[k]


# Filter by transcript expression

## normalizing data (CPM)
input_norm = edger::cpm(GSE125583_Symbol)

## calculating normalized global mean
global_mean <- sum(rowSums(input_norm))/(dim(input_norm)[1] * dim(input_norm)[2])

## removing transcripts with total expression lower than global mean
k <- rowSums(input_norm) >= global_mean
GSE125583_Filtered_Transcripts <- GSE125583_Symbol[k,]


# Filter by sample expression

## removing samples with total expression lower than 10M
k <- colSums(GSE125583_Filtered_Transcripts) >= 10^7
GSE125583_Filtered_Samples <- GSE125583_Filtered_Transcripts[, k]


# Save output data

## this will be repeated later for NC_GSE125583 and All_GSE125583
write.table(GSE125583_Filtered_Samples,
            "~/counts/filtered_data/AD_GSE125583_filtered_samples.txt",
            sep = "\t", col.names = T, row.names = T, quote = F, append = F)
