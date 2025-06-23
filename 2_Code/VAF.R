# Load libraries

library(VariantAnnotation)
library(ggplot2)
library(dplyr)
library(ChIPpeakAnno)
library(ggplot2)
library(RColorBrewer)
library(tibble)
library(reshape2)
library(readxl)
library(grid)
library(tidyverse)

# Set filepath

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma_scWGS/3_Output/Burkitt_scWGS/Figures_supp/")

# Set date

date <- format(Sys.Date(), "%Y%m%d")

# Load functions and color palettes
mycols_paired <- brewer.pal(12,"Paired")
mycols_dark2 <- brewer.pal(8, "Dark2")
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO_Ageline_checks/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Load metadata

below_curve_df <-  read.csv("../../../1_Input/Burkitt_scWGS/below_curve_samples.csv")
bad_baf_df <-  read.csv("../../../1_Input/Burkitt_scWGS/bad_baf_samples.csv")
input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Lymphoma_scWGS/1_Input/Burkitt_scWGS/Sample_overview.xlsx')

# Label input_df with residual info

input_df <- input_df %>%
  mutate(Below_curve = if_else(Sample_name %in% below_curve_df$Sample_name, "Yes", "No"))

# Load vcf files of samples

vcf_files_P3G6_b1 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/P3G6/batch1/snvs/batch1",
                                pattern = "*/*snvs.ptato.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_P3G6_b1 <- vcf_files_P3G6_b1[str_detect(vcf_files_P3G6_b1, "batch1_P3G6GPDABC26|batch1_P3G6GPDABC27")]
vcf_files_P3G6_b2 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/P3G6/batch2/snvs/batch2",
                                pattern = "*/*snvs.ptato.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_P3G6_b2 <- vcf_files_P3G6_b2[!str_detect(vcf_files_P3G6_b2, "batch2_PB11197-BLASC-BCELLP1O3|batch2_PB11197-BLASC-BCELLP1P3")]
vcf_files_P3G6_ptav2 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/P3G6/PTAv2/",
                                   pattern = "*/*snvs.ptatoV2.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_PRN4_b1 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PRN4/batch1/snvs/batch1",
                                pattern = "*/*snvs.ptato.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_PRN4_b2 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PRN4/batch2/snvs/batch2",
                                pattern = "*/*snvs.ptato.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_PRN4_b3 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PRN4/batch3/snvs/batch3",
                                pattern = "*/*snvs.ptato.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_P856_ptav2 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/P856/PTAv2",
                                   pattern = "*/*snvs.ptatoV2.filtered.vcf$", full.names = TRUE, recursive = TRUE )

vcf_files <- c(vcf_files_P3G6_b1, vcf_files_P3G6_b2, vcf_files_P3G6_ptav2, vcf_files_PRN4_b1, vcf_files_PRN4_b2, vcf_files_PRN4_b3, vcf_files_P856_ptav2)
vcf_files <- vcf_files[!str_detect(vcf_files, "PRN4GBDLBC72")] #bulk sample

vcf_list <- lapply(vcf_files, function(f) {
  readVcf(f, genome = "hg38")
})

# Extract VAF per sample function

vaf_df <- do.call(rbind, lapply(seq_along(vcf_list), function(i) {
  vcf <- vcf_list[[i]]
  sample_name <- colnames(vcf)
  
  # Check if VAF field exists and is not empty
  if ("VAF" %in% names(geno(vcf))) {
    vafs <- geno(vcf)$VAF[, 1]  # Extract VAFs for the single sample
    data.frame(Sample_name = sample_name, VAF = as.numeric(vafs))
  } else {
    message(paste("Skipping", sample_name, "- no VAF field"))
    NULL
  }
}))

# Add patient-specific IDs

vaf_df <- left_join(vaf_df, input_df[, c("Sample_name", "Novogene_ID", "BAF", "Below_curve")], 
                    by = c("Sample_name" = "Sample_name"))

# Calculate quantiles and median

median_df <- vaf_df %>%
  group_by(Sample_name) %>%
  summarise(
    Novogene_ID = dplyr::first(Novogene_ID),
    Median = quantile(VAF, 0.5, na.rm = TRUE),
    BAF = dplyr::first(BAF),
    Below_curve = dplyr::first(Below_curve)
  )

# Assign fill color based on residual cutoff

vaf_df$ResidualColor <- ifelse(vaf_df$Below_curve == "Yes",  "lightblue", "lightgrey")
vaf_df$BafColor <- ifelse(vaf_df$BAF == "Bad", "lightgrey",
                             ifelse(vaf_df$BAF == "Intermediate", "lightblue", "deepskyblue4"))

unique_samples <- unique(vaf_df$Novogene_ID)

# Loop over each Patient

for (sample_id in unique_samples) {
  
  sample_data <- vaf_df %>%
    filter(Novogene_ID == sample_id)
  
  sample_median <- median_df %>%
    filter(Novogene_ID == sample_id)
  
  novogene_label <- unique(sample_data$Novogene_ID)
  novogene_label <- ifelse(is.na(novogene_label), sample_id, novogene_label)
  
  ordered_samples <- sample_median %>%
    arrange(desc(Median)) %>%
    pull(Sample_name)
  
  sample_data$Sample_name <- factor(sample_data$Sample_name, levels = ordered_samples)
  
  sample_median$Sample_num <- match(sample_median$Sample_name, levels(sample_data$Sample_name))
  
  p1 <- ggplot(sample_data, aes(x = Sample_name, y = VAF)) +
    geom_violin(aes(fill = ResidualColor), trim = TRUE, color = "black") +
    geom_segment(data = sample_median,
                 aes(x = Sample_num - 0.4, xend = Sample_num + 0.3,
                     y = Median, yend = Median),
                 inherit.aes = FALSE, color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
    scale_fill_identity() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("VAF Distribution for", novogene_label),
         x = "Sample",
         y = "Variant Allele Frequency") +
    coord_cartesian(ylim = c(0, 1))  # safer than using ylim()
  
  ggsave(filename = paste0("VAF_distribution_with_residual_info_", novogene_label, "_", date, ".pdf"),
         plot = p1, width = 5, height = 4)
  
  p2 <- ggplot(sample_data, aes(x = Sample_name, y = VAF)) +
    geom_violin(aes(fill = BafColor), trim = TRUE, color = "black") +
    geom_segment(data = sample_median,
                 aes(x = Sample_num - 0.4, xend = Sample_num + 0.3,
                     y = Median, yend = Median),
                 inherit.aes = FALSE, color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0.4, linetype = "dashed", color = "red") +
    scale_fill_identity() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = paste("VAF Distribution for", novogene_label),
         x = "Sample",
         y = "Variant Allele Frequency") +
    coord_cartesian(ylim = c(0, 1))  # safer than using ylim()
  
  ggsave(filename = paste0("VAF_distribution_with_BAF_info_", novogene_label, "_", date, ".pdf"),
         plot = p2, width = 5, height = 4)
}

# List samples with no VAF=1 variants as these will be probably doublets

vaf1_df <- vaf_df %>%
  group_by(Sample_name) %>%
  summarise(VAF1_count = sum(abs(VAF - 1) < 1e-6)) %>%
  arrange(desc(VAF1_count))

no_vaf1 <- as.character(vaf1_df$Sample_name[vaf1_df$VAF1_count == 0])

no_vaf1_df <- vaf1_df %>%
  filter((Sample_name %in% no_vaf1))

# Samples with low median VAF

median_df <- median_df %>%
  mutate(VAF_low = if_else(Median < 0.4, "Yes", "No"))

low_vaf_df <- median_df %>%
  filter(VAF_low == "Yes")

low_vaf_df <- low_vaf_df %>%
  filter(!(BAF == "Bad" | Below_curve == "Yes"))

# Export samples that did not pass initial QC AND VAF cutoff of 0.4

blacklist_df <- median_df %>%
  filter(BAF == "Bad" | VAF_low == "Yes" | Below_curve == "Yes")

write.csv(blacklist_df, file = "../../../1_Input/Burkitt_scWGS/blacklist_samples.csv", row.names = F)

# Percentage removed because of low quality (poor BAF plot + poor callable loci/mean coverage) + low VAF

perc_removed <- (length(blacklist_df$Sample_name))/length(vcf_files)*100
print(perc_removed)


