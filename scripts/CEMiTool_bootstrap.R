#!/usr/bin/env Rscript
#'
#' This is the R Script to perform filtering step for co-expression analysis.
#'
#' Corresponding authors: Prof. Dr. Gilderlanio Santana de Araujo (gilderlanio [at] gmail.com)
#'                        Ms.C Arthur Ribeiro dos Santos (arthurrdsantos [at] outlook.com)

# Load libraries
library(CEMiTool)
library(dplyr)


# Load AD group data
counts <- read.delim2("~/counts/filtered_data/AD_GSE125583_filtered_samples.txt",
                                             header = T, row.names = 1,
                                             sep = "\t", stringsAsFactors = F)


# Define parameters
## number of repetitions
B = 100L

## number of modules from the original model
nmodules <- 3L

## beta parameter from the original model
beta <- 14L


# Create the initial data frame for number of genes per module
mod_genes <- data.frame(matrix(NA, ncol = nmodules, nrow = B))
colnames(mod_genes) <- c(paste0("M", 1:nmodules))


# Execute loop
for(i in 1:B){
  ## sample counts data
  sample_input <- as.data.frame(t(sample_n(as.data.frame(t(counts)),
                                           214 * 0.7 )))
  
  ## Run cemitool with set_beta for each sampled input
  cem <- cemitool(sample_input, set_beta = beta,
                  apply_vst = TRUE, verbose = FALSE)
  
  ## check number of modules for sampled data
  ### nmodules() also return "Not.Correlated" as a module, subtract result by 1
  nmodules_sample <- (nmodules(cem) - 1L)
  
  ## add columns to mod_genes if nmodules_sample > nmodules
  if(nmodules_sample > nmodules){
    
    ### temporarily store mod_genes
    temp <- mod_genes
    
    ### recreate mod_genes, using a higher number o modules
    mod_genes <- data.frame(matrix(NA, ncol = nmodules_sample, nrow = B))
    
    ### rename columns
    colnames(mod_genes) <- c(paste0("M", 1:nmodules_sample))
    
    ### recover old values
    mod_genes[,1:nmodules] <- temp
    
    ### update nmodules
    nmodules <- nmodules_sample
  }
  
  ## get number of genes per module
  mod_gene_num <- mod_gene_num(cem)
  
  ## add information to mod_genes
  mod_genes[i, 1:nmodules_sample] <- mod_gene_num$num_genes[1:nmodules_sample]
}

# save results
write.table(mod_genes, "~/results/cemitool/Stability/ModulesByBootstrap.txt",
            sep = "\t", col.names = T, row.names = F, quote = F, append = F)