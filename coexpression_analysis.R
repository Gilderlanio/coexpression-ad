#!/usr/bin/env Rscript
#'
#' This is the R Script to perform filtering step for co-expression analysis.
#'
#' Corresponding authors: Prof. Dr. Gilderlanio Santana de Araujo (gilderlanio [at] gmail.com)
#'                        Ms.C Arthur Ribeiro dos Santos (arthurrdsantos [at] outlook.com)
# Load packages
library(CEMiTool)

# Load default pathways
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
paths <- read_gmt(gmt_fname)


# Load default interactions
int_network <- read.delim(system.file("extdata", "interactions.tsv",
                                      package = "CEMiTool"))


# Load data
AD_GSE125583_filtered_samples <- read.delim2("~/counts/filtered_data/AD_GSE125583_filtered_samples.txt",
                                             header = T, row.names = 1,
                                             sep = "\t", stringsAsFactors = F)

NC_GSE125583_filtered_samples <- read.delim2("~/counts/filtered_data/NC_GSE125583_filtered_samples.txt",
                                             header = T, row.names = 1,
                                             sep = "\t", stringsAsFactors = F)

All_GSE125583_filtered_samples <- read.delim2("~/counts/filtered_data/All_GSE125583_filtered_samples.txt",
                                             header = T, row.names = 1,
                                             sep = "\t", stringsAsFactors = F)


# Select group for the process

## this will be repeated later for NC_GSE125583 and All_GSE125583
input_data <- AD_GSE125583_filtered_samples

group_name <- "AD"

# Perform co-expression analysis

cem <- cemitool(input_data, gmt = paths, interactions= int_network,
                apply_vst = TRUE, verbose = TRUE)


# Write files 
generate_report(cem, directory = paste0("~/results/cemitool/", group_name,
                                        "/Report"), force = TRUE)

write_files(cem, directory = paste0("~/results/cemitool/", group_name,
                                    "/Tables"), force = TRUE)

save_plots(cem, "all", directory = paste0("~/results/cemitool/", group_name,
                                          "./Plots"), force = TRUE)


# Repeat process

## change input_data and group_name to match group data
