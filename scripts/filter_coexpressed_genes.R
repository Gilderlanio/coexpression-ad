#!/usr/bin/env Rscript
#'
#' This is the R Script to perform filtering step for co-expression analysis.
#'
#' Corresponding authors: Prof. Dr. Gilderlanio Santana de Araujo (gilderlanio [at] gmail.com)
#'                        Ms.C Arthur Ribeiro dos Santos (arthurrdsantos [at] outlook.com)


# Load libraries
library(edgeR)

# Load data
AD_selected_genes <- read.delim2("~/results/cemitool/AD/Tables/module.tsv",
                                 header = T, sep = "\t", stringsAsFactors = F)

NC_selected_genes <- read.delim2("~/results/cemitool/NC/Tables/module.tsv",
                                 header = T, sep = "\t", stringsAsFactors = F)

All_GSE125583_filtered_samples <- read.delim2("~/counts/filtered_data/All_GSE125583_filtered_samples.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)


# Normalize NC+AD group data by cpm before filtering coexpressed genes
All_GSE125583_cpm <- cpm(All_GSE125583_filtered_samples)


# Filter "Not.Correlated" genes from AD and NC module data
AD_coexpressed_genes <- AD_selected_genes[(AD_selected_genes$modules != "Not.Correlated"), 1]

NC_coexpressed_genes <- NC_selected_genes[(NC_selected_genes$modules != "Not.Correlated"), 1]


# Get unique coexpressed genes
coexpressed_genes <- unique(c(AD_coexpressed_genes, NC_coexpressed_genes))


# Filter normalized NC+AD group by unique coexpressed genes
k <- rownames(All_GSE125583_cpm) %in% unique_coexpressed_genes

All_GSE125583_coexpressed_genes <- All_GSE125583_cpm[k,]


# save output data
write.table(All_GSE125583_coexpressed_genes,
            "~/counts/filtered_data/All_GSE125583_coexpressed_genes.txt",
            sep = "\t", col.names = T, row.names = T, quote = F, append = F)