################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at callable loci and mean coverage of single cell WGS samples (QC step 1)
# Author: Alexander Steemers
# Date: July 2025
################################################################################

# Load libraries

library(ggplot2)
library(tidyverse)
library(readxl)
library(patchwork)
library(RColorBrewer)

# Set filepath

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/")

# Set date

date <- format(Sys.Date(), "%Y%m%d")

# Load plotting functions and color palette

source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO_Ageline_checks/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')
qc_colors <- c('grey', '#54BFB7', '#0A9086')

# Load metadata of single cell samples

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') #dataframe
input_df_sub <- input_df[,c('Sample_name','Mean_coverage','Callable_fraction','Lymphoma_type', "ResolveDNA_version", "CNV", "BAF")] #sub-selection
input_df_sc <- input_df_sub[input_df_sub$ResolveDNA_version %in% c("v1", "v2", "v2.0"), ] #single-cells only
model_df <- input_df_sc[!is.na(input_df_sc$Callable_fraction) & !is.na(input_df_sc$Mean_coverage),] #make sure no N/A in Callable loci and Mean coverage
model_df$Mean_coverage <- as.numeric(model_df$Mean_coverage) #make this numeric

# Decide what cutoff to use for callable loci fraction

ggplot(model_df, aes(x = Mean_coverage, y = Callable_fraction)) +
  geom_point() +
  ylim(c(0,1)) +
  theme_CHemALL() +
  labs(title = "Identification of low quality samples")

ggplot() +
  geom_histogram(data = model_df, aes(x = Callable_fraction),
                 binwidth = 0.06)+ 
  geom_vline(xintercept = 0.45) + 
  theme_CHemALL() +
  scale_fill_manual(values = qc_colors[c(1,3)]) +
  ggtitle('Callable loci distribution')

cf_cutoff <- 0.45 # decided based on histogram above
below_cf_cutoff <- subset(model_df, Callable_fraction < cf_cutoff) 
model_df$CallableLociQuality <- 'Pass'
model_df[model_df$Sample_name %in% below_cf_cutoff$Sample_name,]$CallableLociQuality <- 'Fail'

# Plot histogram

ggplot() +
  geom_histogram(data = model_df, aes(x = Callable_fraction, fill = CallableLociQuality),
                 binwidth = 0.06)+ 
  geom_vline(xintercept = 0.45) + 
  theme_CHemALL() +
  scale_fill_manual(values = qc_colors[c(1,3)]) +
  ggtitle('Callable loci distribution')
ggsave(paste0("Figures/callable_loci_histogram_allsamples_", date, ".pdf"), width = 5, height = 3)

# Remove samples with lower than 40% callable genomes

model_df_sub <- subset(model_df, Callable_fraction >= 0.45) # remove samples with lower than 40% callable genomes

# Samples which are initially filtered out due to low callable fraction (dataframe to export later)

low_call_frac_df <- subset(model_df, Callable_fraction < 0.45)

# Now plot logistic curve for mean coverage vs callable loci fraction
# Manual predicted values, define start values here

L_start <- max(model_df_sub$Callable_fraction)
k_start <- 0.1
x0_start <- mean(model_df_sub$Mean_coverage) 

# Fit logistic model

nls_logistic <- nls(Callable_fraction ~ L / (1 + exp(-k * (Mean_coverage - x0))),
                    data = model_df_sub,
                    start = list(L = L_start, k = k_start, x0 = x0_start))

# Add predictions

model_df_sub$predicted_logistic <- predict(nls_logistic)

# Plot logistic curve

ggplot(model_df_sub, aes(x = Mean_coverage, y = Callable_fraction)) +
  geom_point(color = "#54BFB7") +
  scale_y_continuous(limits = c(0.40, 1)) +
  geom_line(aes(y = predicted_logistic), color = "#000000", size = 1) +
  labs(title = "Logistic Fit to CallableFraction vs MeanCoverage")

# Calculate residuals to identify outliers

model_df_sub$residual <- model_df_sub$Callable_fraction - model_df_sub$predicted_logistic
cutoff <- -0.05 # decided based on histogram below
below_curve <- subset(model_df_sub, residual < cutoff)
model_df_sub$MappingQuality <- 'Pass'
model_df_sub[model_df_sub$Sample_name %in% below_curve$Sample_name,]$MappingQuality <- 'Fail'

# Plot residuals in histogram and set cutoff for low quality samples

ggplot() +
  geom_histogram(data = model_df_sub, aes(x = residual, fill = MappingQuality),
                 binwidth = 0.02)+ 
  geom_vline(xintercept = cutoff) + 
  theme_CHemALL() +
  scale_fill_manual(values = qc_colors[c(1,3)]) +
  ggtitle('Residual from predicted callable genome fraction')
ggsave(paste0("Figures/residual_histogram_allsamples_", date, ".pdf"), width = 5, height = 3)

# Plot original relationship, but now with mappping quality filter

ggplot(model_df_sub, aes(x = Mean_coverage, y = Callable_fraction)) +
  theme_CHemALL() +
  geom_line(aes(y = predicted_logistic), color = qc_colors[2], linewidth = 1) +
  geom_point(aes(color = MappingQuality), size = 2) +
  scale_color_manual(values = qc_colors[c(1,3)]) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7))+
  labs(title = "Identification of low quality samples")
ggsave(paste0("Figures/correlation_plot_allsamples_", date, ".pdf"), width = 5, height = 3)

# Print outlier samples

below_curve$Sample_name

# Export samples that did not pass initial QC 

write.csv(below_curve, file = "Data/below_curve_samples.csv", row.names = F)
write.csv(low_call_frac_df , file = "Data/low_callable_loci.csv")

# Side analysis

# Plot category CNV and BAFplot vs residual

model_df_plot <- model_df_sub[ !(is.na(model_df_sub$BAF)),]
model_df_plot$CNV <- factor(model_df_plot$CNV, levels = c('to do','Bad','Intermediate','Good'))
model_df_plot$BAF <- factor(model_df_plot$BAF, levels = c('to do','Bad','Intermediate','Good'))

p1 <- ggplot(data = model_df_plot[model_df_plot$CNV != 'to do',], aes(x = CNV, y = residual, fill = CNV)) +
  geom_violin() + scale_fill_manual(values = qc_colors) + theme_CHemALL() + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  geom_hline(yintercept = cutoff, linetype = "solid", color = "black") +
  ggtitle('Residual from predicted callable genome fraction')

p2 <- ggplot(data = model_df_plot[model_df_plot$BAF != 'to do',], aes(x = BAF, y = residual, fill = BAF)) +
  geom_violin() + scale_fill_manual(values = qc_colors) + theme_CHemALL() + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  geom_hline(yintercept = cutoff, linetype = "solid", color = "black") +
  ggtitle('Residual from predicted callable genome fraction')

p1 + p2
ggsave(paste0("Figures/residual_density_perCNV_BAF_group_", date, ".pdf"), width = 7, height = 3)

p3 <- ggplot(data = model_df_plot[model_df_plot$BAF != 'to do',], aes(x = BAF, y = Callable_fraction, fill = BAF)) +
  geom_violin() + scale_fill_manual(values = qc_colors) + theme_CHemALL() + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.7)

p4 <- ggplot(data = model_df_plot[model_df_plot$BAF != 'to do',], aes(x = BAF, y = Mean_coverage, fill = BAF)) +
  geom_violin() + scale_fill_manual(values = qc_colors) + theme_CHemALL() + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.7)

p3 + p4

# Intersect between bad BAF and low callable samples

bad_baf_df <- subset(input_df, BAF == "Bad")

intersect(bad_baf_df$Sample_name, below_curve$Sample_name)
setdiff(bad_baf_df$Sample_name, below_curve$Sample_name)
low_q_samples <- union(bad_baf_df$Sample_name, below_curve$Sample_name)
low_q_samples_df <- subset(model_df_sub, Sample_name %in% low_q_samples)
other_filters <- unique(c(below_curve$Sample_name, low_call_frac_df$Sample_name))
unique_to_bad_baf <- setdiff(bad_baf_df$Sample_name, other_filters)





# I want to check whether the samples that branch off earlier than the major clone are low-quality

samples_to_label <- c("PRN4GPDLBC15", "PRN4GPDLBC17", "PB08410-BLBM-BCELLP2G8", "PVA9GTDABC55", "")

ggplot(model_df_sub, aes(x = Mean_coverage, y = Callable_fraction)) +
  geom_line(aes(y = predicted_logistic), color = qc_colors[2], linewidth = 1) +
  geom_point(aes(color = MappingQuality), size = 2) +
  geom_text(data = filter(model_df_sub, Sample_name %in% samples_to_label),
            aes(label = Sample_name), vjust = -1, size = 3) +
  scale_color_manual(values = qc_colors[c(1,3)]) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7)) +
  theme_CHemALL() +
  labs(title = "Identification of low quality samples")
