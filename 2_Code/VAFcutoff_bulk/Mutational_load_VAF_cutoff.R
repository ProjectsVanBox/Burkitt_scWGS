################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to save mutational load and decide on vaf cutoffs for bulk WGS samples 
# Author: Alexander Steemers
# Date: July 2025
################################################################################

# Load libraries

library(VariantAnnotation)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)

# Set working directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/VAF")

# Load functions and plotting functions

source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Samples 

Patients <- c("P3G6", "PRN4", "P856", "PJBU", "PIA9", "PVA9")
Patients <- c("P3G6")

# File path template

vcf_base <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/"
#vcf_suffix <- "/HaplotypeCaller/vcf/germline/gatk4haplotypecaller/"
vcf_suffix <- "/HaplotypeCaller/vcf/germline/somatic_filtering/SMuRF/"

# Empty list to collect data

# Steps: Take only PASS variants --> only SNVs --> only autosomal --> VAF >= 0.15

all_patients_data <- list()

get_ref <- function(x) if (!is.null(x) && length(x) >= 1) as.numeric(x[1]) else NA
get_alt <- function(x) if (!is.null(x) && length(x) >= 2) as.numeric(x[2]) else NA

for (Patient in Patients) {
  
  #vcf_file_unzipped <- paste0(vcf_base, Patient, vcf_suffix, Patient, ".vep.vcf")
  #vcf_file_gzipped <- paste0(vcf_base, Patient, vcf_suffix, Patient, ".vep.vcf.gz")
  
  vcf_file_unzipped <- paste0(vcf_base, Patient, vcf_suffix, Patient, "_bulk.vep.SMuRF.filtered.sorted.VAF03.vcf")
  vcf_file_gzipped <- paste0(vcf_base, Patient, vcf_suffix, Patient, "_bulk.vep.SMuRF.filtered.sorted.VAF03.vcf.gz")
  
  if (file.exists(vcf_file_unzipped)) {
    vcf_file <- vcf_file_unzipped
  } else if (file.exists(vcf_file_gzipped)) {
    vcf_file <- vcf_file_gzipped
  } else {
    cat("No VCF file found for", Patient, "\n")
    next
  }
  
  vcf <- readVcf(vcf_file, "hg38")
  
  # Keep only PASS variants
  vcf <- vcf[as.vector(unlist(fixed(vcf)$FILTER)) == "PASS"]
  
  # Keep only SNVs (remove indels)
  ref_allele <- as.character(ref(vcf))
  alt_list <- alt(vcf)
  
  is_snv <- mapply(function(ref, alt) {
    nchar(ref) == 1 && all(nchar(as.character(alt)) == 1)
  }, ref_allele, alt_list)
  
  #vcf <- vcf[is_snv]
  
  num_pass <- length(vcf)
  cat("Patient:", Patient, "- PASS variants:", num_pass, "\n")
  
  # Autosomes only
  vcf <- vcf[seqnames(rowRanges(vcf)) %in% as.character(1:22)]
  
  # REF / ALT depths
  AD <- geno(vcf)$AD
  ref_depth <- matrix(sapply(AD, get_ref), nrow = dim(AD)[1], ncol = dim(AD)[2])
  alt_depth <- matrix(sapply(AD, get_alt), nrow = dim(AD)[1], ncol = dim(AD)[2])
  
  vaf <- alt_depth / (ref_depth + alt_depth)
  vaf[vaf < 0.15 | is.na(vaf)] <- NA
  
  # Variant names
  variant_ids <- paste0(seqnames(rowRanges(vcf)), ":", start(rowRanges(vcf)))
  
  # Prepare dataframe
  vaf_df <- as.data.frame(vaf)
  colnames(vaf_df) <- colnames(vcf)
  vaf_df$Variant <- variant_ids
  
  # Long format
  vaf_long <- pivot_longer(vaf_df, cols = -Variant, names_to = "Sample", values_to = "VAF") %>%
    filter(!is.na(VAF)) %>%
    mutate(Donor = Patient)
  
  all_patients_data[[Patient]] <- vaf_long
  
  cat("Done with Patient:", Patient, "\n")
}

# Combine all patients

combined_df <- do.call(rbind, all_patients_data)
PRN4_P1B11_df <- combined_df[combined_df$Sample == "PB08410-BLLN-BCELLP1B11", ]
P3G6_P1B4_df <- combined_df[combined_df$Sample == "PB11197-BLASC-BCELLP1B4", ]
P3G6_msc_df <- combined_df[combined_df$Sample == "PB11197-BLBM-MSCBULK", ]

# Define bulk samples

bulk_WGS <- c("PB11197-BLASC-BCELLBULK", "PRN4GBDLBC72", "PB14458-BLPL-BCELLBULK","PB14458-BLBM-BCELLBULK", "PIA9GBDABC78", "PVA9GBDABC78","PJBUGBDABC82")

# Filter out MSCs

combined_df <- combined_df %>%
  filter(!grepl("MS", Sample))

# Separate bulk and single-cell

bulk_df <- combined_df %>% filter(Sample %in% bulk_WGS)
single_cell_df <- combined_df %>% filter(!Sample %in% bulk_WGS)

# Violin plots with cutoffs for bulk 30X samples

for (donor in Patients) {
  
  donor_df <- bulk_df %>% filter(Donor == donor)
  
  # Decide intercept
  if ("PRN4GBDLBC72" %in% donor_df$Sample) {
    intercept <- 0.25
  } else {
    intercept <- 0.3
  }
  
  p <- ggplot(donor_df, aes(x = Sample, y = VAF, fill = Sample)) +
    geom_violin(trim = FALSE, scale = "width", fill = "steelblue") +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white", color = "black") +
    theme_bw() +
    geom_hline(yintercept = intercept, colour = "red") +
    scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
  theme(axis.text.x = element_text(),
          strip.text = element_text(size = 10),
          strip.background = element_rect(fill = "grey90", color = NA)) +
    labs(title = paste0("Patient ", donor, " - Bulk sample(s)"),
         y = "VAF", x = "Sample") +
    guides(fill = "none")
  
  print(p)
  
  ggsave(paste0("Figures/VAF_distribution_bulk_", donor, ".pdf"), plot = p, width = 7, height = 4)
}

# Apply bulk-specific VAF cutoffs (determined by violin plots)

bulk_df <- bulk_df %>%
  rowwise() %>%
  mutate(VAF = ifelse(
    Sample == "PRN4GBDLBC72" & VAF < 0.25, NA,
    ifelse(Sample != "PRN4GBDLBC72" & VAF < 0.3, NA, VAF)
  )) %>%
  filter(!is.na(VAF)) %>%
  ungroup()

# Calculate purity of bulks (determine which are uncontaminated and which are contaminated)

bulk_purity <- bulk_df %>%
  group_by(Sample) %>%
  summarise(median_VAF = median(VAF, na.rm = TRUE)) %>%
  arrange(desc(median_VAF))

bulk_purity$Purity <- bulk_purity$median_VAF*2 # PRN4 LN contaminated --> adjust VAF cutoff
print(bulk_purity)

# Violin plots with cutoffs for single cells 15X samples

single_cell_df <- single_cell_df %>%
  mutate(Sample = factor(Sample),
         Donor = factor(Donor))

unique_donors <- unique(single_cell_df$Donor)

for (donor in unique_donors) {
  
  donor_df <- single_cell_df %>% 
    filter(Donor == donor)
  
  # Order samples by median VAF (high to low)
  sample_order <- donor_df %>%
    group_by(Sample) %>%
    summarise(median_VAF = median(VAF, na.rm = TRUE)) %>%
    arrange(desc(median_VAF)) %>%
    pull(Sample)
  
  donor_df <- donor_df %>%
    mutate(Sample = factor(Sample, levels = sample_order))
  
  p <- ggplot(donor_df, aes(x = Sample, y = VAF)) +
    geom_violin(fill = "steelblue", color = "black", trim = FALSE, scale = "width") +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white", color = "black") +
    theme_bw() +
    geom_hline(yintercept = 0.15, color = "red") +
    scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
    theme(
      axis.text.x = element_text(size = 4,angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    ) +
    labs(
      title = paste0("Patient ", donor, " - Single-cell samples"),
      y = "VAF",
      x = "Sample"
    )
  
  print(p)
  
  ggsave(paste0("Figures/VAF_distribution_single_cell_", donor, ".pdf"), plot = p, width = 7, height = 4)
  
}

# Save RDS files

saveRDS(bulk_df, "../MutLoad/Data/PASS_autosomal_VAFcutoff_bulk_samples.RDS")
saveRDS(single_cell_df, "../MutLoad/Data/PASS_autosomal_VAFcutoff_single_cell_samples.RDS")

# Optional read RDS here

bulk_df <- readRDS(file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Data/PASS_autosomal_VAFcutoff_bulk_samples.RDS")
single_cell_df <- readRDS(file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Data/PASS_autosomal_VAFcutoff_single_cell_samples.RDS")
          
# Remove below curve samples and determine VAF cutoff using MAD outlier detection

below_curve_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv")
low_call_frac_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_frac_samples.csv")

single_cell_df_filtered <- single_cell_df %>%
  filter(!Sample %in% below_curve_df$Sample_name,
         !Sample %in% low_call_frac_df)

# Parameters

num_bins <- 10
epsilon <- 1e-6  # to avoid zero-probability issues

# Decide here if using filtered data (coverage quality) or unfiltered data (all samples)

# 1. Bin VAFs into histogram for each sample
binned_df <- single_cell_df_filtered %>%
  mutate(bin = cut(VAF, breaks = seq(0, 1, length.out = num_bins + 1), include.lowest = TRUE)) %>%
  group_by(Sample, bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(prob = (count + epsilon) / sum(count + epsilon)) %>%  # normalize
  ungroup()

# 2. Pivot to wide format: one row per sample
wide_probs <- binned_df %>%
  dplyr::select(Sample, bin, prob) %>%
  pivot_wider(names_from = bin, values_from = prob, values_fill = list(prob = epsilon)) %>%
  column_to_rownames("Sample")

# 3. Compute reference distribution (e.g., median across samples)
ref_dist <- apply(wide_probs, 2, median)

# 4. Compute TVD for each sample
tvd <- function(p, q) {
  0.5 * sum(abs(p - q))
}
tvd_values <- apply(wide_probs, 1, function(p) tvd(p, ref_dist))

# 5. Output: samples ranked by TVD
tvd_df <- data.frame(Sample = names(tvd_values), TVD = tvd_values) %>%
  arrange(desc(TVD))

# 6. Plot TVD scores
mad_val <- mad(tvd_df$TVD)
median_val <- median(tvd_df$TVD)

tvd_df <- tvd_df %>%
  mutate(Flagged = TVD > (median_val + 2.5 * mad_val)) # https://www.sciencedirect.com/science/article/pii/S0022103113000668?via%3Dihub

ggplot(tvd_df, aes(x = TVD, y = reorder(Sample, TVD), 
                   fill = Flagged)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Total Variation Distance of VAF Distributions",
       x = "TVD", y = "Sample_name") +
  scale_fill_manual(values = c('#54BFB7','grey')) + theme_CHemALL() +
  theme(text =  element_text(size =  7, color = 'black'),
        axis.text = element_text(angle = 90, size = 5, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5))

# Annotate the other plot with this
plot_df3b <- merge(single_cell_df, tvd_df)

# Remove flagged with higher than median VAF
median(plot_df3b$VAF)
median_df <- plot_df3b %>% group_by(Sample) %>% summarise(med = median(VAF))

plot_df3b[plot_df3b$Sample %in% median_df[median_df$med > median(plot_df3b$VAF),]$Sample, 'Flagged'] <- FALSE

# Rename that column
plot_df3b$VAFfilter <- 'Pass'
plot_df3b[plot_df3b$Flagged,]$VAFfilter <- 'Fail'

unique_donors <- unique(plot_df3b$Donor)

for (donor in unique_donors) {
  
  donor_df <- plot_df3b %>% 
    filter(Donor == donor)
  
  # Calculate median VAF per sample
  sample_medians <- donor_df %>% 
    group_by(Sample) %>% 
    summarise(median_vaf = median(VAF, na.rm = TRUE)) %>% 
    arrange(desc(median_vaf))
  
  # Reorder Sample factor based on median VAF
  donor_df$Sample <- factor(donor_df$Sample, levels = sample_medians$Sample)
  
  # Plot
  p <- ggplot(data = donor_df,
              aes(x = Sample,
                  y = VAF,
                  fill = VAFfilter)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.8) +
    ggtitle(paste0('VAF distribution per single-cell sample - ', donor)) +
    scale_fill_manual(values = c('Fail' = 'grey', 'Pass' = '#54BFB7')) +
    theme_CHemALL() +
    ggTextAxisRotate() +
    theme(text = element_text(size = 7, color = 'black'),
          axis.text = element_text(size = 5, colour = "black"))
  
  print(p)
}

fail_samples <- plot_df3b %>%
  filter(VAFfilter == "Fail") %>%
  pull(Sample) %>%
  unique() %>%
  as.character()

blacklist_samples <- unique(c(filtered_out, fail_samples, below_curve_df$Sample_name))

P856_blacklist <- blacklist_samples[grepl("^P856|^PB14458", blacklist_samples)]
P3G6_blacklist <- blacklist_samples[grepl("^P3G6|^PB11197", blacklist_samples)]
PRN4_blacklist <- blacklist_samples[grepl("^PRN4|^PB08410", blacklist_samples)]
PIA9_blacklist <- blacklist_samples[grepl("^PIA9", blacklist_samples)]
PVA9_blacklist <- blacklist_samples[grepl("^PVA9", blacklist_samples)]
PJBU_blacklist <- blacklist_samples[grepl("^PJBU", blacklist_samples)]

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx')
colnames(input_df)
input_df_sc <- input_df[input_df$ResolveDNA_version %in% c("v1", "v2", "v2.0"), ]

perc_removed <- (length(blacklist_samples))/length(input_df_sc$Sample_name)*100
print(perc_removed)

abs_num_kept <- length(input_df_sc$Sample_name) - length(blacklist_samples)

input_df_sub <- input_df_sc[!(input_df_sc$Sample_name %in% blacklist_samples), ]

input_df_sub %>%
  group_by(Novogene_ID) %>%
  summarise(Number_of_sc = n())
