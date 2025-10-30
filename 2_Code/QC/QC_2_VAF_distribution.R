################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at VAF distribution of single cell WGS samples (QC step 2)
# Author: Alexander Steemers
# Date: July 2025
################################################################################

# Laod libraries and functions

library(tibble)   
library(dplyr)   
library(tidyr)   
library(ggplot2) 
library(reshape2)
library(VariantAnnotation)
library(readxl)
library(ggpubr)
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Set filepath

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC")

# Load autosomal PASS variants above 0.15 VAF

SBSs_PASS <- readRDS(file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Data/autosomal_PASS_variants_VAF015.RDS")
INDELs_PASS <- readRDS(file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Data/autosomal_INDEL_PASS_variants_VAF015.RDS")

# load metadata info

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx')
colnames(input_df)
input_df_sub <- input_df[,c('Sample_name','Mean_coverage' ,'Novogene_ID','Callable_fraction','Lymphoma_type', "ResolveDNA_version", "CNV", "BAF")]
input_df_sc <- input_df_sub[input_df_sub$ResolveDNA_version %in% c("v1", "v2", "v2.0"), ]
below_curve_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv")
low_call_frac_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv")
bad_baf_df <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv")

# Generate initial blacklist sample (based on QC step 1) and remove remaining bad baf samples and filter those out

blacklist_samples <- unique(c(below_curve_df$Sample, low_call_frac_df$Sample_name, bad_baf_df$Sample_name)) # filtering on bad baf as well here

input_df_sc_filtered <- input_df_sc %>%
  filter(!Sample_name %in% blacklist_samples)

# Plot TVD for all single cell mutations

# NB: variants underwent the following selection: SMuRF --> PTATO --> autosomal --> VAF >= 0.15

vafCheck <- SBSs_PASS[names(SBSs_PASS) %in% input_df_sc_filtered$Sample_name]

# Convert granges to list of dataframes

empty_df <- list()
for (Sample in names(vafCheck)){
  vafCheck[[Sample]]$samplename <- Sample
  
  print(head(data.frame(vafCheck[[Sample]])))
  
  empty_df[[Sample]] <- data.frame(vafCheck[[Sample]])
}

plot_df <- bind_rows(empty_df, .id = "column_label")
plot_df1 <- merge(plot_df, input_df_sc_filtered[,c('Sample_name','Novogene_ID')], by.x = 'samplename',by.y = 'Sample_name')

# Determine cutoff using MAD outlier detection

# Parameters

num_bins <- 10
epsilon <- 1e-6  # to avoid zero-probability issues

# Decide here if using filtered data (coverage quality) or unfiltered data (all samples) --> I used filtered data

input_df_vafs <- plot_df1 # choose here
input_df_vafs <- input_df_vafs %>%
  mutate(VAF = as.numeric(VAF))

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

# Plot TVD of VAF distributions

ggplot(tvd_df, aes(x = TVD, y = reorder(samplename, TVD), 
                   fill = Flagged)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Total Variation Distance of VAF Distributions",
       x = "samplename", y = "TVD") +
  scale_fill_manual(values = c('#54BFB7','grey')) + theme_CHemALL() +
  theme(text =  element_text(size =  7, color = 'black'),
        axis.text = element_text(size = 2, colour = "black")) +
  ggTextAxisRotate()
ggsave('Figures/TVD_ranked_flagged_3_mad.pdf', width = 10, height = 3)

# Annotate the other plot with this

plot_df3b <- merge(input_df_vafs, tvd_df)

# Remove flagged with higher than median VAF

median(plot_df3b$VAF)
median_df <- plot_df3b %>% group_by(samplename) %>% summarise(med = median(VAF))

plot_df3b[plot_df3b$samplename %in% median_df[median_df$med > median(plot_df3b$VAF),]$samplename, 'Flagged'] <- FALSE

# Rename that column

plot_df3b$VAFfilter <- 'Pass'
plot_df3b[plot_df3b$Flagged,]$VAFfilter <- 'Fail'

# Split the dataframe into a list by Novogene_ID

plot_list <- split(plot_df3b, plot_df3b$Novogene_ID)

# Create one ggplot per Novogene_ID

plot_list <- lapply(names(plot_list), function(id) {
  ggplot(data = plot_list[[id]],
         aes(x = samplename,
             y = VAF,
             fill = VAFfilter)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.8) +
    ggtitle(paste('VAF distribution -', id)) +
    scale_fill_manual(values = c('Fail' = 'grey', 'Pass' = '#54BFB7')) +
  theme_CHemALL() +
    ggTextAxisRotate() +
    theme(text = element_text(size = 7, color = 'black'),
          axis.text = element_text(size = 4, colour = "black"))
})
names(plot_list) <- names(split(plot_df3b, plot_df3b$Novogene_ID))  # optional, to keep names

plot_list[[1]]  
plot_list[[2]]  
plot_list[[3]]  
plot_list[[4]] 
plot_list[[5]] 
plot_list[[6]]  

for (i in seq_along(plot_list)) {
  plot_name <- names(plot_list)[i]
  
  ggsave(filename = paste0("Figures/VAF_violinplot_TVD_flagged_", plot_name, ".pdf"),
         plot = plot_list[[i]],
         width = 6, height = 4, units = "in")
}

# Export these samples that failed VAF QC step 2

fail_df_pta <- unique(plot_df3b[plot_df3b$VAFfilter == 'Fail', c('samplename','Novogene_ID')])
write_csv(fail_df_pta, '../QC/Data/PTA_samples_failVAFcheck.txt')

other_filters <- unique(c(below_curve$Sample_name, low_call_frac_df$Sample_name, fail_df_pta$Sample_name))
unique_to_bad_baf <- setdiff(bad_baf_df$Sample_name, other_filters)

# Side analysis

filtered_samples <- unique(c(low_call_frac_df$Sample_name, below_curve_df$Sample_name, bad_baf_df$Sample_name, fail_df_pta$samplename))  # samples that didn't pass QC

length(input_df_sc$Sample_name) - length(filtered_samples) # number of cells left after filtering steps 1 and 2 
((length(input_df_sc$Sample_name) - length(filtered_samples))/length(input_df_sc$Sample_name)) *100 # percentage of cells left


pb11197_p3g6 <- filtered_samples[grepl("^P3G6|^PB11197", filtered_samples)]
pb08410_prn4 <- filtered_samples[grepl("^PB08410|^PRN4", filtered_samples)]
pb14458_p856 <- filtered_samples[grepl("^PB14458|^P856", filtered_samples)]
pia9        <- filtered_samples[grepl("^PIA9", filtered_samples)]
pva9        <- filtered_samples[grepl("^PVA9", filtered_samples)]
pjbu        <- filtered_samples[grepl("^PJBU", filtered_samples)]

group_counts <- data.frame(
  Group = c("PB11197 / P3G6", "PB08410 / PRN4", "PB14458 / P856", "PIA9", "PVA9", "PJBU"),
  Count = c(length(pb11197_p3g6),
            length(pb08410_prn4),
            length(pb14458_p856),
            length(pia9),
            length(pva9),
            length(pjbu))
)

print(group_counts) # how many cells were removed because of QC steps

sample_names <- as.character(input_df_sc$Sample_name)

# Get all sample names
all_samples <- unique(input_df_sc$Sample_name)

# Get samples that passed QC
kept_samples <- setdiff(all_samples, filtered_samples)

# Count kept samples per group
pb11197_p3g6_kept <- kept_samples[grepl("^P3G6|^PB11197", kept_samples)]
pb08410_prn4_kept <- kept_samples[grepl("^PB08410|^PRN4", kept_samples)]
pb14458_p856_kept <- kept_samples[grepl("^PB14458|^P856", kept_samples)]
pia9_kept         <- kept_samples[grepl("^PIA9", kept_samples)]
pva9_kept         <- kept_samples[grepl("^PVA9", kept_samples)]
pjbu_kept         <- kept_samples[grepl("^PJBU", kept_samples)]

# Create a dataframe of counts
group_counts_kept <- data.frame(
  Group = c("PB11197 / P3G6", "PB08410 / PRN4", "PB14458 / P856", "PIA9", "PVA9", "PJBU"),
  Kept = c(length(pb11197_p3g6_kept),
           length(pb08410_prn4_kept),
           length(pb14458_p856_kept),
           length(pia9_kept),
           length(pva9_kept),
           length(pjbu_kept))
)

print(group_counts_kept)

median(group_counts_kept$Kept) 
mean(group_counts_kept$Kept) 


# Get some metadata now from the cells that passed all QC
input_df_QC_pass <- input_df_sc %>% filter(!Sample_name %in% filtered_samples)
input_df_QC_pass$Mean_coverage <- as.numeric(input_df_QC_pass$Mean_coverage) #make this numeric

# Get mean genomes per patient (+ range)
summary(input_df_QC_pass$Mean_coverage)

# Get total autosomal SNVs and INDELS
filtered_SBSs_PASS <- SBSs_PASS[!names(SBSs_PASS) %in% filtered_samples]
print(sum(sapply(filtered_SBSs_PASS, length)))

filtered_INDELs_PASS <- INDELs_PASS[!names(INDELs_PASS) %in% filtered_samples]
print(sum(sapply(filtered_INDELs_PASS, length)))

