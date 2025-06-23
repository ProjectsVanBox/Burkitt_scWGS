################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at callable loci, coverage and BAF quality of single cell WGS samples 
# Written by: Alexander Steemers
# Date: June 2025
# Modified: 
################################################################################

# Load libraries

library(ggplot2)
library(tidyverse)
library(readxl)
library(patchwork)

# Set filepath

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/")

# Set date

date <- format(Sys.Date(), "%Y%m%d")

# Load plotting functions and color palette

source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO_Ageline_checks/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')
qc_colors <- c('grey', '#54BFB7', '#0A9086')

# Load metadata of samples

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx')
colnames(input_df)
input_df_sub <- input_df[,c('Sample_name','Mean_coverage','Callable_fraction','Lymphoma_type', "ResolveDNA_version", "CNV", "BAF")]
input_df_sc <- input_df_sub[input_df_sub$ResolveDNA_version %in% c("v1", "v2", "v2.0"), ]

# Fit logaritmic model to Callable vs Mean_Coverage to identify outliers (i.e. likely lower quality)

model_df <- input_df_sc[!is.na(input_df_sc$Callable_fraction) & !is.na(input_df_sc$Mean_coverage),]

# Manual predicted values, define start values here

L_start <- max(model_df$Callable_fraction)
k_start <- 0.5
x0_start <- mean(model_df$Mean_coverage) 

# Fit logistic model

nls_logistic <- nls(Callable_fraction ~ L / (1 + exp(-k * (Mean_coverage - x0))),
                    data = model_df,
                    start = list(L = L_start, k = k_start, x0 = x0_start))

# Add predictions

model_df$predicted_logistic <- predict(nls_logistic)

# Plot

ggplot(model_df, aes(x = Mean_coverage, y = Callable_fraction)) +
  geom_point(color = "#54BFB7") +
  geom_line(aes(y = predicted_logistic), color = "#000000", size = 1) +
  labs(title = "Logistic Fit to CallableFraction vs MeanCoverage")

# Calculate residuals to identify outliers

model_df$residual <- model_df$Callable_fraction - model_df$predicted_logistic
cutoff <- -0.05
below_curve <- subset(model_df, residual < cutoff)
model_df$MappingQuality <- 'Pass'
model_df[model_df$Sample_name %in% below_curve$Sample_name,]$MappingQuality <- 'Fail'

#### plot residuals in histogram and set cutoff for low quality samples

ggplot() +
  geom_histogram(data = model_df, aes(x = residual, fill = MappingQuality),
                 binwidth = 0.1) + 
  geom_vline(xintercept = cutoff) + 
  theme_CHemALL() +
  scale_fill_manual(values = qc_colors[c(1,3)]) +
  ggtitle('Residual from predicted callable genome fraction')
ggsave(paste0("QC/residual_histogram_allsamples_", date, ".pdf"), width = 5, height = 3)

# Plot original relationship

ggplot(model_df, aes(x = Mean_coverage, y = Callable_fraction)) +
  theme_CHemALL() +
  geom_line(aes(y = predicted_logistic), color = qc_colors[2], linewidth = 1) +
  geom_point(aes(shape = Lymphoma_type, color = MappingQuality)) +
  #geom_text(aes(label = Sample_ID), hjust = 0, vjust = 1.2, size = 2) +
  scale_color_manual(values = qc_colors[c(1,3)]) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7))+
  labs(title = "Identification of low quality samples")
ggsave(paste0("QC/correlation_plot_allsamples_", date, ".pdf"), width = 5, height = 3)

# Print outlier samples

below_curve$Sample_name

# Plot category CNV and BAFplot vs residual

model_df_plot <- model_df[ !(is.na(model_df$BAF)),]
model_df_plot$CNV <- factor(model_df_plot$CNV, levels = c('to do','Bad','Intermediate','Good'))
model_df_plot$BAF <- factor(model_df_plot$BAF, levels = c('to do','Bad','Intermediate','Good'))

p1 <- ggplot(data = model_df_plot[model_df_plot$CNV != 'to do',], aes(x = CNV, y = residual, fill = CNV)) +
  geom_violin() + scale_fill_manual(values = qc_colors) + theme_CHemALL() + 
  geom_hline(yintercept = cutoff, linetype = "solid", color = "black") +
  ggtitle('Residual from predicted callable genome fraction')

p2 <- ggplot(data = model_df_plot[model_df_plot$BAF != 'to do',], aes(x = BAF, y = residual, fill = BAF)) +
  geom_violin() + scale_fill_manual(values = qc_colors) + theme_CHemALL() + 
  geom_hline(yintercept = cutoff, linetype = "solid", color = "black") +
  ggtitle('Residual from predicted callable genome fraction')

p1 + p2
ggsave(paste0("QC/residual_density_perCNV_BAF_group_", date, ".pdf"), width = 7, height = 3)

# Intersect between bad BAF and low callable samples

bad_baf_df <- subset(input_df, BAF == "Bad")
intersect(bad_baf_df$Sample_name, below_curve$Sample_name)
setdiff(bad_baf_df$Sample_name, bad_baf_df$Sample_name)
low_q_samples <- union(bad_baf_df$Sample_name, bad_baf_df$Sample_name)
low_q_samples_df <- subset(model_df, Sample_name %in% low_q_samples)

# Percentage removed because of low quality (poor BAF plot + poor callable loci/mean coverage)

perc_removed <- (length(low_q_samples))/length(model_df$Sample_name)*100
print(perc_removed)

# Export samples that did not pass initial QC 

write.csv(below_curve, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/below_curve_samples.csv", row.names = F)
write.csv(bad_baf_df, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/bad_baf_samples.csv", row.names = F)
