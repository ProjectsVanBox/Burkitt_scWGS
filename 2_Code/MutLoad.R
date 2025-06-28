################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to get filtered SNVs (total + autosomal) from single cell WGS samples 
# Author: Alexander Steemers
# Date: June 2025
################################################################################

# Load libraries

library(reshape2)
library(ggplot2)
library(tidyverse)
library(VariantAnnotation)
library(readxl)
library(BSgenome)
library(GenomicRanges)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)

# Set working directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad")

# Load functions and plotting functions

source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Load metadata

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx')
input_df_sc <- input_df[input_df$ResolveDNA_version %in% c("v1", "v2", "v2.0"), ]
input_df_sub <- input_df_sc[!is.na(input_df_sc$Callable_fraction) & !is.na(input_df_sc$Mean_coverage),]

# Define PTATO directory

ptato_dir <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO"

# Which patients to check

folders_to_check <- c("P3G6", "PRN4", "P856")

# MinimalVAF cut-off (0.15 is used for PTA data so will stick to this cutoff for now)

MinimalVAF <- 0.20

# Function to locate the right vcf file because some are V2 filtered and others are not

find_filtered_vcf <- function(sample_id, ptato_dir, subdirs) {
  patterns <- c(
    sprintf("%s.*\\.snvs\\.ptatoV2\\.filtered\\.vcf(\\.gz)?$", sample_id),  # preferred
    sprintf("%s.*\\.snvs\\.ptato\\.filtered\\.vcf(\\.gz)?$",   sample_id)   # fallback
  )
  for (pat in patterns) {
    for (d in subdirs) {
      hit <- list.files(file.path(ptato_dir, d),
                        pattern     = pat,
                        recursive   = TRUE,
                        full.names  = TRUE,
                        ignore.case = TRUE)
      if (length(hit)) return(hit[1])
    }
  }
  NA_character_
}

# Loop over all samples and make a list 

SBSs_filtered  <- list()   # every variant, with FILTER label

for (Sample in input_df_sub$Sample_name) {
  message("→ processing ", Sample)
  
  vcf_path <- find_filtered_vcf(Sample, ptato_dir, folders_to_check)
  if (is.na(vcf_path)) {
    warning("No filtered PTATO VCF found for ", Sample)
    next
  }
  
  print(vcf_path) # to check if the right VCF file was used
  
  vcf <- readVcf(vcf_path)
  
  # Read VAF
  if (!"VAF" %in% names(geno(vcf)))
    stop("`geno(vcf)$VAF` not present in ", basename(vcf_path))
  
  vaf <- geno(vcf)$VAF
  if (length(dim(vaf)) == 2L)       
    vaf <- vaf[, 1, drop = TRUE]
  
  # Add FILTER logic
  rr <- rowRanges(vcf)
  mcols(rr)$Chromosome <- as.character(seqnames(rr))
  mcols(rr)$VAF    <- vaf
  mcols(rr)$FILTER <- ifelse(vaf > MinimalVAF, "PASS", "FAIL_VAF")
  
  SBSs_filtered[[Sample]]  <- rr
}

# Save both autosomal and sex chromosome SNVs

saveRDS(SBSs_filtered, file = "Data/filtered_variants_ResolveDNA_VAF020.RDS")

# Possibility to load from here

#SBSs_filtered <- readRDS(file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Data/filtered_variants_ResolveDNA_VAF015.RDS")

# Split vaf > minimalvaf and vaf =< minimal vaf

SBSs_PASS_filtered_total <- lapply(SBSs_filtered, function(x) x[which(x$FILTER =="PASS"),])
SBSs_FAIL_VAF_filtered_total <- lapply(SBSs_filtered, function(x) x[which(x$FILTER =="FAIL_VAF"),])

# Save R objects 

saveRDS(SBSs_PASS_filtered_total, file = "Data/total_filtered_PASS_variants_PTA_VAF020.RDS")
saveRDS(SBSs_FAIL_VAF_filtered_total, file = "Data/total_filtered_FAIL_VAF_variants_PTA_VAF020.RDS")

# Only autosomal SNVs

auto_chrs <- as.character(1:22)

SBSs_filtered_autosomal <- lapply(SBSs_filtered, function(gr) {
    keep <- seqnames(gr) %in% auto_chrs
    gr2  <- gr[keep]                        
    keepSeqlevels(gr2, auto_chrs, pruning.mode = "coarse")
  }
)

# Save R object because parts above take very long to load

saveRDS(SBSs_filtered_autosomal, file = "Data/autosomal_filtered_variants_ResolveDNA_VAF020.RDS")

# Possibility to load from here

#SBSs_filtered_autosomal <- readRDS(file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Data/autosomal_filtered_variants_ResolveDNA_VAF015.RDS")

# Remove all variants that are uncallable and/or with a PTAprobs < PTAprobsCutoff

SBSs_PASS_filtered_autosomal <- lapply(SBSs_filtered_autosomal, function(x) x[which(x$FILTER =="PASS"),])
SBSs_FAIL_VAF_filtered_autosomal <- lapply(SBSs_filtered_autosomal, function(x) x[which(x$FILTER =="FAIL_VAF"),])

# Save R objects 

saveRDS(SBSs_PASS_filtered_autosomal, file = "Data/autosomal_PASS_variants_PTA_VAF020.RDS")
saveRDS(SBSs_FAIL_VAF_filtered_autosomal, file = "Data/autosomal_FAIL_VAF_variants_PTA_VAF020.RDS")

# Get lengths for each SBS label per sample

get_lengths <- function(lst) vapply(lst, length, integer(1))

counts <- tibble(
  Sample    = names(SBSs_PASS_filtered_autosomal),                     
  PASS      = get_lengths(SBSs_PASS_filtered_autosomal),
  FAIL      = get_lengths(SBSs_FAIL_VAF_filtered_autosomal)
) |>
  pivot_longer(-Sample, names_to = "Category", values_to = "Count")

# Plot stacked barplot to see how much is filtered out

ggplot(counts, aes(x = Sample, y = Count, fill = Category)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = c(PASS = "forestgreen", FAIL = "firebrick")) +
  labs(x = "Sample", y = "Number of autosomal SNVs", fill = "Filter status") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(paste0("Figures/Mut_load_filtering_stacked_barplot_VAF020.pdf"))
