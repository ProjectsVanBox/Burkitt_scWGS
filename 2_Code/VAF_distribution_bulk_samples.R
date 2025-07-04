################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at VAF distributions of single cell WGS samples
# Author: Alexander Steemers
# Date: June 2025
################################################################################

# Load libraries

library(VariantAnnotation)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tibble)
library(readxl)
library(tidyverse)
library(purrr)

# Set filepath

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/VAF")

# Set date

date <- format(Sys.Date(), "%Y%m%d")

# Load functions and color palettes

mycols_paired <- brewer.pal(12,"Paired")
mycols_dark2 <- brewer.pal(8, "Dark2")
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO_Ageline_checks/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Load metadata

diagnostic_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_overview.csv')

# Load bulk WGS Burkitt VCF files from diagnostics (SMuRF filtered as these are PASS + VAF > 0.3?)

bulkWGS_diagnostic_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples",
                                           pattern = "*/*/*filtered.sorted.vcf.gz$", full.names = TRUE, recursive = TRUE)
bulkWGS_diagnostic_vcf_files_sub <- bulkWGS_diagnostic_vcf_files[!str_detect(bulkWGS_diagnostic_vcf_files, "PMCID211AAO.vep")] # filter out duplicate

# Get IDs of tumour columns for each diagnostic sample

tumour_ids <- unique(diagnostic_df$Tumor)

# Extract VAFs of variants for each sample

vaf_tbl_list <- map(bulkWGS_diagnostic_vcf_files_sub, function(vcf_path) {
  
  vcf <- readVcf(vcf_path)                                  
  
  sample_cols <- intersect(colnames(vcf), tumour_ids) 
  if (length(sample_cols) == 0) {
    message("Skipping ", vcf_path, " (no matching tumour ID)")
    return(NULL)
  }
  
  if (!"VAF" %in% names(geno(vcf))) {           
    message("Skipping ", vcf_path, " (no VAF field)")
    return(NULL)
  }
  
  rr <- rowRanges(vcf)
  ref_allele <- as.character(mcols(rr)$REF)
  alt_allele <- as.character(unlist(mcols(rr)$ALT))
  var_id <- paste0(seqnames(rr), ":", start(rr), "_", ref_allele, ">", alt_allele)
  
  vaf_mat <- geno(vcf)$VAF
  
  out <- tibble(variant = var_id)
  for (sid in sample_cols) {
    out[[sid]] <- as.numeric(vaf_mat[, sid])
  }
  out
})

# Combine all tibbles into one wide data.frame

vaf_df <- reduce(vaf_tbl_list, full_join, by = "variant")

# Order from highest median VAF to lowest median VAF

vaf_long <- vaf_df %>%
  pivot_longer(-variant,
               names_to  = "Sample",
               values_to = "VAF") %>%
  filter(!is.na(VAF)) %>%
  group_by(Sample) %>%
  mutate(median_VAF = median(VAF, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Sample = reorder(Sample, -median_VAF))

# Make violin plot

pdf("Figures/VAF_distribution_diagnostic_tumor_samples.pdf", width = 10, height = 6)
ggplot(vaf_long, aes(x = Sample, y = VAF)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.10, outlier.size = 0.5) +  # optional inset box
  geom_hline(yintercept = 0.25,linetype = "dashed", colour = "red") +
  geom_hline(yintercept = 0.20,linetype = "dashed", colour = "blue") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(size = 14, face = "bold")) +
  ylab("Variant Allele Frequency") +
  xlab("Tumour sample") +
  ggtitle("VAF distributions across diagnostic tumour samples")
dev.off()

# Inspect range of VAF values

range(vaf_long$VAF)

# Get median VAFs and then decipher the tumour purity

median_tbl <- vaf_long %>%                                   
  group_by(Sample) %>%
  summarise(
    median_VAF        = median(VAF, na.rm = TRUE),
    purity_estimate   = 2 * median_VAF,                      
    .groups = "drop"
  ) %>%
  mutate(
    recommended_cutoff = case_when(
      purity_estimate >  0.70                     ~ 0.25,    
      purity_estimate >= 0.60 & purity_estimate <= 0.70 ~ 0.20,    
      TRUE                                           ~ NA_real_    
    )
  )

# Get a summary

high_purity   <- median_tbl %>% filter(recommended_cutoff == 0.25) %>% pull(Sample)
mid_purity    <- median_tbl %>% filter(recommended_cutoff == 0.20) %>% pull(Sample)

cat("Cut-off 0.25 (>70% purity):",
    if (length(high_purity)) paste(high_purity, collapse = ", ") else "none", "\n")
cat("Cut-off 0.20 (60–70% purity):",
    if (length(mid_purity))  paste(mid_purity,  collapse = ", ") else "none", "\n")

# Load bulk WGS Burkitt VCF files

vcf_file_PRN4_bulk_lymph_node <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PRN4/batch1/snvs/batch1/PB08410_PRN4.vep_PRN4GBDLBC72/PB08410_PRN4.vep_PRN4GBDLBC72.snvs.ptato.filtered.vcf.gz"
vcf_file_P3G6_bulk_ascites <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB11197_new/intermediate/short_variants/somatic_vcfs/PB11197/231018_HMFreg2086_PB11197.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_PB11197-BLASC-BCELLBULK.SMuRF.filtered.vcf.gz"
vcf_file_P856_bulk_pleura <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB14458/intermediate/short_variants/somatic_vcfs/PB14458/231123_HMFreg2101_PB14458.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_PB14458-BLPL-BCELLBULK.SMuRF.filtered.vcf.gz"
vcf_file_P856_bulk_bone_marrow <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB14458/intermediate/short_variants/somatic_vcfs/PB14458/231123_HMFreg2101_PB14458.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_PB14458-BLBM-BCELLBULK.SMuRF.filtered.vcf.gz"
bulkWGS_vcf_files <- c(vcf_file_PRN4_bulk_lymph_node, vcf_file_P3G6_bulk_ascites, vcf_file_P856_bulk_pleura, vcf_file_P856_bulk_bone_marrow)

# Make function to extract variant info

variant_key <- function(rr) {
  paste0(
    seqnames(rr), ":", start(rr), "_",
    as.character(mcols(rr)$REF), ">",       
    as.character(unlist(mcols(rr)$ALT))        
  )
}

# Build one tibble per VCF, then full-join

vaf_tbl_list2 <- map(bulkWGS_vcf_files, function(vcf_path) {
  
  vcf <- readVcf(vcf_path)
  
  if (!"VAF" %in% names(geno(vcf))) {
    message("Skipping ", vcf_path, " (no VAF field)")
    return(NULL)
  }
  
  rr       <- rowRanges(vcf)
  var_id   <- variant_key(rr)
  vaf_mat  <- geno(vcf)$VAF         
  samples  <- colnames(vcf)         
  
  out <- tibble(variant = var_id)
  for (sid in samples) {
    out[[sid]] <- as.numeric(vaf_mat[, sid])
  }
  out
})

# Merge all per-sample tables into one wide data frame

vaf_df2 <- reduce(vaf_tbl_list2, full_join, by = "variant")

vaf_long2 <- vaf_df2 %>%
  pivot_longer(-variant,
               names_to  = "Sample",
               values_to = "VAF") %>%
  filter(!is.na(VAF))          # drop missing calls

# Make violin plot 

pdf("Figures/VAF_distribution_bulk_sorted_tumor_samples.pdf",
    width = 10, height = 6)
ggplot(vaf_long2, aes(x = Sample, y = VAF)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.10, outlier.size = 0.5) +  
  geom_hline(yintercept = 0.25,linetype = "dashed", colour = "red") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1)) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(size = 14, face = "bold")) +
  ylab("Variant Allele Frequency") +
  xlab("Tumour sample") +
  ggtitle("VAF distributions across bulk tumour samples")
dev.off()

# Inspect range of VAF values

range(vaf_long2$VAF)
