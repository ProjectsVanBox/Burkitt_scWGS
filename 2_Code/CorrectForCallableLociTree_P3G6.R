################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to pick CellPhy tree and correct for callable loci 
# Author: Alexander Steemers
# Date: July 2025
################################################################################

# Set directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/")

# Load libraries

library(tidyverse)
library(reshape2)
library(ggtree)
library(treeio)
library(stringi)
library(dplyr)
library(cellPhyWrapperPlotting)  
library(phangorn)
library(ggplot2)
library(patchwork)
library(readxl)

# Define functions to make callable loci-corrected trees

# Function to calculate the product of (1 - sensitivity) for each sample in a node
calculate_product_sensitivity <- function(samples, df_sensitivity) {
  # Split the 'samples' string by '|'
  sample_list <- strsplit(samples, "\\|")[[1]]
  
  # Retrieve sensitivities for these samples
  sensitivities <- df_sensitivity$sensitivity[df_sensitivity$Sample_name %in% sample_list]
  
  # Calculate the product of (1 - sensitivity)
  1 - prod(1 - sensitivities)
}

# function to take a tree and correct the branch lengths based on the product sensitivity
## by calculating the product of (1 - sensitivity) for each sample in a node
# requires samples (=samples column from output of cellphywrapper object: tree@data$samples .
# requires a CALLABLE df (=containing columns SAMPLES and CALLABLE (in nr of bases))

correct_branches <- function(tree, callable_df){
  # get sensitivity df by taking fraction of total possible callable loci
  max_callable <- 2745186691
  callable_df$sensitivity <- callable_df$Callable_Loci / max_callable
  sensitivity_df <- callable_df[c('Sample_name','sensitivity')]
  
  print(sensitivity_df)
  
  # calculate the sensitivity per node/branch using the calculate_product_sensitivity function
  tree@data$product_sensitivity <- sapply(tree@data$samples, calculate_product_sensitivity, sensitivity_df)
  
  # correct the branch lengths
  tree@data$corr_branch_lengths <- tree@data$branch_length / tree@data$product_sensitivity
  
  print(tree@data$branch_length)
  print(tree@data$corr_branch_lengths)
  
  tree@data$corr_branch_lengths[is.na(tree@data$corr_branch_lengths)] <- 0 # for merged germline sample
  
  # store also in the phylo object, required for downstream analyses based on the phylo object
  tree@phylo$edge.length <- tree@data[order(tree@data$node),]$corr_branch_lengths
  
  # store in branch_length column for cellphyplotting use
  tree@data$branch_length <- tree@data$corr_branch_lengths
  
  print(tree@phylo$edge.length)
  
  return(tree)
}

# Import callable loci and SBS PASS for donor

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') #dataframe
input_df_P3G6 <- input_df[input_df$Novogene_ID == "P3G6", ]
exclude_list <- readLines("~/hpc/pmc_vanboxtel/projects/Burkitt/2_Code/CellPhyWrapper/P3G6_v2/exclude_bulks_and_blacklisted_samples.txt")
input_df_P3G6_filtered <- input_df_P3G6[!(input_df_P3G6$Sample_name %in% exclude_list) & !grepl("MSC", input_df_P3G6$Sample_name), ]
input_df_P3G6_filtered$Callable_Loci <- as.numeric(input_df_P3G6_filtered$Callable_Loci)

SBSs_autosomal_PASS <- readRDS(file = "../MutLoad/Data/autosomal_PASS_variants_VAF015.RDS")
filtered_SBS_list <- SBSs_autosomal_PASS[names(SBSs_autosomal_PASS) %in% input_df_P3G6_filtered$Sample_name]

mut_load_df <- data.frame(
  Sample_name        = names(filtered_SBS_list),   
  Number_of_mutations = lengths(filtered_SBS_list), 
  row.names           = NULL,
  stringsAsFactors    = FALSE
)

# Print each CellPhy iteration

for (i in c(1,2,3,4,5,6,7,8,9)) {
  # Construct the path to the current tree file
  tree_path <- paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/P3G6/Filtered_samples/CPW_0", i, "/TreeObject0.", i, ".RDS")
  
  # Read the tree object
  tree <- readRDS(tree_path)
  
  # Correct the branches
  corrected_tree <- correct_branches(tree, input_df_P3G6_filtered)
  
  # Plot before and after
  plot_before <- plot_gg_tree(tree)
  plot_after <- plot_gg_tree(corrected_tree)
  
  # Combine and title
  combined_plot <- (plot_before / plot_after) + plot_annotation(title = paste("TreeObject0.", i, sep = ""))
  
  # Display
  print(combined_plot)
  
  # Optionally compare trees
  #comparePhylo(tree@phylo, corrected_tree@phylo, force.rooted = TRUE, plot = TRUE, use.edge.length = FALSE)
}

normal_samples <- c(
  "PB11197-BLASC-BCELLP2C4",
  "PB11197-BLASC-BCELLP2D4",
  "PB11197-BLASC-BCELLP2E4",
  "PB11197-BLASC-BCELLP2B4",
  "PB11197-BLASC-BCELLP2F4",
  "P3G6GPDABC31"
)

hg38_autosomal_nonN_genome_size <- 2745186691
mut_load_df$SNV_LOAD_NORM <- mut_load_df$Number_of_mutations/input_df_P3G6_filtered$Callable_Loci*hg38_autosomal_nonN_genome_size

# Calculate averages
snv_avg_summary <- mut_load_df %>%
  mutate(group = if_else(Sample_name %in% normal_samples, "normal", "tumor")) %>%
  group_by(group) %>%
  summarise(avg_snv_load = mean(SNV_LOAD_NORM, na.rm = TRUE))

snv_avg_summary

# Conclusion: Percentage 20% seems to be the best

tree_chosen <- paste0("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Markus/P3G6/CPW_0", 4, "/TreeObject0.", 4, ".RDS")

# Read the tree object
tree <- readRDS(tree_chosen)

# Correct the branches
corrected_tree <- correct_branches(tree, callable_df_P3G6)
corrected_tree@data$branch_length <- round(corrected_tree@data$branch_length)

# Plot final tree

pdf("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Markus/P3G6/Corrected_absolute_tree_20perc.pdf", width = 10, height = 6)  # adjust size as needed

plot_gg_tree(corrected_tree, add_tip_label = F, add_title = "P3G6 | Ascites | 13.7Y", add_branch_length = T)

dev.off()
