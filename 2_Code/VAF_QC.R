################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to check for VAF QC in single cell WGS samples
# Author: Alexander Steemers
# Date: June 2025
################################################################################

# Load libraries

library(tibble)   
library(dplyr)   
library(tidyr)   
library(ggplot2) 
library(reshape2)
library(VariantAnnotation)
library(readxl)
library(ggpubr)
library(readr)

# Set working directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/VAF")

# Load functions and plotting functions

source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Load preprocessed PTATO output (filtered for VAF, PTAprob)

SBSs_PASS <- readRDS(file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Data/autosomal_PASS_variants_PTA_VAF015.RDS")
#SBSs_FAIL <- readRDS(file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Data/autosomal_FAIL_variants_PTA_VAF015.RDS")

# Load metadata info and filter

input_df <- read_excel(path = '~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx')
input_df_sc <- input_df[input_df$ResolveDNA_version %in% c("v1", "v2", "v2.0"), ]
input_df_sub <- input_df_sc[!is.na(input_df_sc$Callable_fraction) & !is.na(input_df_sc$Mean_coverage),]

# Open blacklist (mapping quality) and filter out

below_curve_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv")
bad_baf_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv")
low_quality_samples <- unique(c(below_curve_df$Sample_name, bad_baf_df$Sample_name))
input_df_filter <- input_df_sub %>% mutate(LowQual = ifelse(Sample_name %in% low_quality_samples, "Yes", "No"))

# Plot now for all FAIL mutations > what is the distribution there? (these are smurf filtered)

#vafCheckfail <- SBSs_FAIL[names(SBSs_PASS) %in% input_df_filter$Sample_name]

# Convert granges to list of dataframes

#empty_df <- list()
#for (Sample in names(vafCheckfail)){
#  vafCheckfail[[Sample]]$samplename <- Sample
#  
#  print(head(data.frame(vafCheckfail[[Sample]])))
#  
#  empty_df[[Sample]] <- data.frame(vafCheckfail[[Sample]])
#}

#plot_dffail <- bind_rows(empty_df, .id = "column_label")

# Add cell type and donor for faceting

#plot_df2 <- merge(plot_dffail, input_df_filter[,c('Novogene_ID','Myc_translocation_IGV', 'Sample_name', 'ResolveDNA_version', 'LowQual')], by.x = 'samplename',by.y = 'Sample_name')
#plot_df2filt <- plot_df2[plot_df2$LowQual != "Yes", ]

#ggplot(data = plot_df2filt,
#       aes(y = VAF, x = Myc_translocation_IGV)) +
#  geom_abline(slope = 0, intercept = 0.4) +
#  geom_violin(draw_quantiles = c(0.5, 0.75, 0.95), alpha = 0.8) +
#  scale_fill_manual(values = c('grey', '#54BFB7', '#0A9086')) +
#  ggtitle('VAF distribution post-PTATO') +
#  theme_CHemALL() +
#  theme(text =  element_text(size =  7, color = 'black'),
#        axis.text = element_text(size = 5, colour = "black"))
# No big difference between myc translocated and non-translocated B cells within the FAIL variants

#ggplot(data = plot_df2,
#       aes(y = VAF, x = LowQual)) +
#  geom_abline(slope = 0, intercept = 0.4) +
#  geom_violin(draw_quantiles = c(0.5, 0.75, 0.95), alpha = 0.8) +
#  scale_fill_manual(values = c('grey', '#54BFB7', '#0A9086')) +
#  ggtitle('VAF distribution post-PTATO') +
#  theme_CHemALL() +
#  theme(text =  element_text(size =  7, color = 'black'),
#        axis.text = element_text(size = 5, colour = "black"))
# No big difference between lowqual and highqual samples within the FAIL variants

#median_quantile <- quantile(plot_df2filt$VAF)[[3]]
#threequart_quantile <- quantile(plot_df2filt$VAF)[[4]]

# Plot now for all PASS mutations > what is the distribution there? (these are smurf filtered)

vafCheck <- SBSs_PASS[names(SBSs_PASS) %in% input_df_filter$Sample_name]

# convert granges to list of dataframes
empty_df <- list()
for (Sample in names(vafCheck)){
  vafCheck[[Sample]]$samplename <- Sample
  
  print(head(data.frame(vafCheck[[Sample]])))
  
  empty_df[[Sample]] <- data.frame(vafCheck[[Sample]])
}

plot_df <- bind_rows(empty_df, .id = "column_label")

# Add cell type and donor for faceting

plot_df3 <- merge(plot_df, input_df_filter[,c('Novogene_ID','Myc_translocation_IGV', 'Sample_name', 'ResolveDNA_version', 'LowQual')], by.x = 'samplename',by.y = 'Sample_name')
plot_df3filt <- plot_df3[plot_df3$LowQual != "Yes", ]

plot_df3filt <- plot_df3filt %>% mutate(
    VAF = vapply(VAF, function(x) as.numeric(x[1]), numeric(1L))
  )

# Plot distribution per ResolveDNA version

ggplot(data = plot_df3filt,
       aes(y = VAF, x = ResolveDNA_version)) +
  geom_abline(slope = 0, intercept = 0.4) +
  geom_violin(draw_quantiles = c(0.5), alpha = 0.8) +
  scale_fill_manual(values = c('grey', '#54BFB7', '#0A9086')) +
  ggtitle('VAF distribution post-PTATO') +
  theme_CHemALL() +
  theme(text =  element_text(size =  7, color = 'black'),
        axis.text = element_text(size = 5, colour = "black"))
# v1 samples seem to have a slightly better distribution

# And per Donor

ggplot(data = plot_df3filt,
       aes(y = VAF, x = samplename, fill = FILTER)) +
  geom_abline(slope = 0, intercept = 0.4) +
  geom_hline(yintercept = 0.15, linetype = "dashed", colour = "red") +
  geom_violin(draw_quantiles = c(0.5), alpha = 0.8) +
  scale_fill_manual(values = c( '#0A9086','#54BFB7','grey')) +
  ggtitle('VAF distribution post-PTATO') +
  theme_CHemALL() +
  ylim(0,1) + 
  ggTextAxisRotate() +
  facet_wrap(~ Novogene_ID, scales = "free", ncol = 3) +
  theme(text =  element_text(size =  7, color = 'black'),
        axis.text = element_text(size = 5, colour = "black"))

# Merge and plot together

#merge_plot_df <- rbind(plot_df2filt, plot_df3filt)

#ggplot(data = merge_plot_df,
#       aes(x = FILTER,
#           y = VAF, fill = FILTER)) +
  #geom_abline(slope = 0, intercept = 0.4, color = '#54BFB7') +
  #geom_abline(slope = 0, intercept = median_quantile, color = 'grey') +
  #geom_abline(slope = 0, intercept = threequart_quantile, color = 'grey') +
#  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.8) +
#  scale_fill_manual(values = c('grey', '#54BFB7', '#0A9086')) +
#  ggtitle('VAF distribution per sample post-PTATO') +
  #facet_wrap(~ CellType + DONOR, scales = "free", ncol = 7) +
#  theme_CHemALL() +
#  ylim(0,1) +
#  ggTextAxisRotate() +
#  theme(text =  element_text(size =  7, color = 'black'),
#        axis.text = element_text(size = 5, colour = "black"))

#ggsave('Figures/PTA_samples_PASS-FAIL_postPTATOVAFdistr.pdf',
#       width = 75, height = 75, units = 'mm')

#ggplot(data = merge_plot_df,
#       aes(x = samplename,
#           y = VAF, fill = FILTER)) +
  #geom_abline(slope = 0, intercept = 0.4, color = '#54BFB7') +
  #geom_abline(slope = 0, intercept = median_quantile, color = 'grey') +
  #geom_abline(slope = 0, intercept = threequart_quantile, color = 'grey') +
#  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.8) +
#  scale_fill_manual(values = c('grey', '#54BFB7', '#0A9086')) +
#  ggtitle('VAF distribution per sample post-PTATO') +
#  facet_wrap(~ samplename + Novogene_ID, scales = "free", ncol = 12) +
#  theme_CHemALL() +
#  ggTextAxisRotate() +
#  theme(text =  element_text(size =  7, color = 'black'),
#        axis.text = element_text(size = 5, colour = "black"))#

#ggsave('Figures/PTA_samples_postPTATOVAFdistr.pdf',
#       width = 300, height = 300, units = 'mm')


# Plot only PASS samples and determine cutoff using MAD outlier detection

# Parameters

num_bins <- 10
epsilon <- 1e-6  # to avoid zero-probability issues

# Decide here if using filtered data (coverage quality) or unfiltered data (all samples)

input_df_vafs <- plot_df3filt # choose here

# 1. Bin VAFs into histogram for each sample
binned_df <- input_df_vafs %>%
  mutate(bin = cut(VAF, breaks = seq(0, 1, length.out = num_bins + 1), include.lowest = TRUE)) %>%
  group_by(samplename, bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(samplename) %>%
  mutate(prob = (count + epsilon) / sum(count + epsilon)) %>%  # normalize
  ungroup()

# 2. Pivot to wide format: one row per sample
wide_probs <- binned_df %>%
  dplyr::select(samplename, bin, prob) %>%
  pivot_wider(names_from = bin, values_from = prob, values_fill = list(prob = epsilon)) %>%
  column_to_rownames("samplename")

# 3. Compute reference distribution (e.g., median across samples)
ref_dist <- apply(wide_probs, 2, median)

# 4. Compute TVD for each sample
tvd <- function(p, q) {
  0.5 * sum(abs(p - q))
}
tvd_values <- apply(wide_probs, 1, function(p) tvd(p, ref_dist))

# 5. Output: samples ranked by TVD
tvd_df <- data.frame(samplename = names(tvd_values), TVD = tvd_values) %>%
  arrange(desc(TVD))

# 6. Plot TVD scores
mad_val <- mad(tvd_df$TVD)
median_val <- median(tvd_df$TVD)

tvd_df <- tvd_df %>%
  mutate(Flagged = TVD > (median_val + 3.5 * mad_val)) # https://www.sciencedirect.com/science/article/pii/S0022103113000668?via%3Dihub

ggplot(tvd_df, aes(x = TVD, y = reorder(samplename, TVD), 
                   fill = Flagged)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Total Variation Distance of VAF Distributions",
       x = "TVD", y = "Sample_name") +
  scale_fill_manual(values = c('#54BFB7','grey')) + theme_CHemALL() +
  theme(text =  element_text(size =  7, color = 'black'),
        axis.text = element_text(size = 5, colour = "black")) +
  ggTextAxisRotate()
ggsave('Figures/PTA_samples_TVD_ranked_flagged_35mad.pdf', width = 10, height = 3)

# Annotate the other plot with this
plot_df3b <- merge(input_df_vafs, tvd_df)

# Remove flagged with higher than median VAF
median(plot_df3b$VAF)
median_df <- plot_df3b %>% group_by(samplename) %>% summarise(med = median(VAF))

plot_df3b[plot_df3b$samplename %in% median_df[median_df$med > median(plot_df3b$VAF),]$samplename, 'Flagged'] <- FALSE

# Rename that column
plot_df3b$VAFfilter <- 'Pass'
plot_df3b[plot_df3b$Flagged,]$VAFfilter <- 'Fail'

ggplot(data = plot_df3b,
       aes(x = samplename,
           y = VAF,
           fill = VAFfilter)) +
  #geom_abline(slope = 0, intercept = 0.4, color = '#54BFB7') +
  #geom_abline(slope = 0, intercept = median_quantile, color = 'grey') +
  #geom_abline(slope = 0, intercept = threequart_quantile, color = 'grey') +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.8) +
  ggtitle('VAF distribution per sample post-PTATO') +
  facet_grid(~ Novogene_ID, scales = "free", space = 'free') +
  scale_fill_manual(values = c('grey','#54BFB7')) + theme_CHemALL() +
  ggTextAxisRotate() +
  theme(text =  element_text(size =  7, color = 'black'),
        axis.text = element_text(size = 5, colour = "black"))
ggsave('Figures/PTA_postPTATOvafs_TVD_flagged_35mad.pdf', width = 12, height = 5)

# Export these failed VAF samples

fail_df_pta <- unique(plot_df3b[plot_df3b$VAFfilter == 'Fail', c('samplename','Novogene_ID')])
write_csv(fail_df_pta, 'Data/PTA_samples_failVAFcheck_35mad.txt')

