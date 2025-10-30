################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to do O/E ratio
# Author: Alexander Steemers
# Date: September 2025
################################################################################

# Load libraries
library(MutationalPatterns)
library(BSgenome)
library(VariantAnnotation)
library(ChIPpeakAnno)
library(ggplot2)
library(NMF)
library(RColorBrewer)
library(nlme)
library(tibble)
library(reshape2)
library(grid)
library(readxl)
library(stringr)
library(tidyr)
library(ggpubr)
library(dplyr)
library(VariantAnnotation)
library(purrr)
library(forcats)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# Load functions and colour palettes

mycols_paired <- brewer.pal(12,"Paired")
mycols_dark2  <- brewer.pal(8, "Dark2")
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Set directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/")

# Load metadata 

input_df       <- read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
diagnostic_df  <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_overview.csv')
low_callable_df<- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv')
below_curve_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df     <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df    <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')

# Make blacklist

blacklist <- unique(c(below_curve_df$Sample_name,
                      low_callable_df$Sample_name,
                      bad_baf_df$Sample_name,
                      fail_vaf_df$samplename))

# Load single cell VCF files (PTATO filtered)

ptato_dir <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO"
folders_to_check <- c("P3G6", "PRN4", "P856", "PIA9", "PVA9", "PJBU")

all_filtered_vcfs <- unlist(lapply(folders_to_check, function(subdir) {
  list.files(
    file.path(ptato_dir, subdir),
    pattern = "snvs.*filtered\\.vcf(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )
}))

# Remove old PTATO file names

all_filtered_vcfs <- all_filtered_vcfs[!grepl("old", all_filtered_vcfs, ignore.case = TRUE)]
  
# Samples re-run with PTAv2 filtering → drop earlier filtered files

samples_to_exclude <- c(
  "PB11197-BLASC-BCELLP1B4","PB11197-BLASC-BCELLP1C4","PB11197-BLASC-BCELLP1I4",
  "PB11197-BLASC-BCELLP1J3","PB11197-BLASC-BCELLP1K4","PB11197-BLASC-BCELLP1L3",
  "PB11197-BLASC-BCELLP1O3","PB11197-BLASC-BCELLP1P3","P3G6GDDABC71",
  "PB14458-BLPL-BCELLP4B3","PB14458-BLPL-BCELLP4B5","PB14458-BLPL-BCELLP4C3",
  "PB14458-BLPL-BCELLP4D3","PB14458-BLPL-BCELLP4D5","PB14458-BLPL-BCELLP4E3",
  "PB14458-BLPL-BCELLP4J3","PB14458-BLPL-BCELLP4K3","PB14458-BLPL-BCELLP4K5",
  "PB14458-BLPL-BCELLP4L3","PB14458-BLPL-BCELLP4L5","PB14458-BLPL-BCELLP4M3",
  "P856GDDUBC32","P856GDDUBC33","P856GDDUBC34","P856GDDUBC40","P856GDDUBC41",
  "P856GDDUBC42","P856GDDUBC43","P856GDDUBC44","P856GDDUBC45",
  "PB14458-BLBM-BCELLP2B3","PB14458-BLBM-BCELLP2B4","PB14458-BLBM-BCELLP2C4",
  "PB14458-BLBM-BCELLP2E4","PB14458-BLBM-BCELLP2F2","PB14458-BLBM-BCELLP2F4",
  "PB14458-BLBM-BCELLP2I2","PB14458-BLBM-BCELLP2L3","PB14458-BLBM-BCELLP2L4",
  "PB14458-BLBM-BCELLP2M4","PB14458-BLBM-BCELLP2N2","PB14458-BLBM-BCELLP2N4",
  "P856GDDBBC46","P856GDDBBC48","P856GDDBBC54","P856GDDBBC57","P856GDDBBC58",
  "P856GDDBBC59","P856GDDBBC60","P856GDDBBC61","P856GDDBBC62","P856GDDBBC63","P856GDDBBC64"
)

pattern_exclude <- paste0("(", paste0(samples_to_exclude, collapse = "|"), ").*\\.snvs\\.ptato\\.filtered\\.vcf\\.gz$")
all_filtered_vcfs <- all_filtered_vcfs[!grepl(pattern_exclude, all_filtered_vcfs, ignore.case = TRUE)]

# use non blacklisted samples vcf!!!!
# Filter out blacklist samples

single_cell_sample_names <- sub(".*\\.vep_([^/\\.]+).*", "\\1", all_filtered_vcfs)
scWGS_vcf_files_sub      <- all_filtered_vcfs[!single_cell_sample_names %in% blacklist]
single_cell_sample_names_sub <- single_cell_sample_names[!single_cell_sample_names %in% blacklist]

# Read in SNVs per genome

grl <- read_vcfs_as_granges(scWGS_vcf_files_sub, single_cell_sample_names_sub, genome = "hg38", type = "snv")

# Define autosomes
autosomes <- paste0("chr", 1:22)

# Keep only autosomes in the GRangesList
grl_autosomes <- endoapply(grl, function(x) x[seqnames(x) %in% autosomes])

snv_counts <- sapply(grl_autosomes, length)

# Turn into data frame

snv_df <- data.frame(
  Sample_name = names(snv_counts),
  SNV_count = as.numeric(snv_counts),
  stringsAsFactors = FALSE
)

# Correct for callable loci

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') #dataframe

# Merge with input_df by Sample
merged_df <- merge(input_df, snv_df, by = "Sample_name", all.x = TRUE)

# Normalize counts by callable fraction
merged_df$SNV_per_callable <- merged_df$SNV_count / merged_df$Callable_fraction

filtered_df <- merged_df[, c("Sample_name", "Age_at_sampling_Y", "SNV_per_callable", "Myc_translocation_IGV", "Biopsy_type" )]

sbs9_status_df <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/Normal_SBS9_status_table.csv", stringsAsFactors =  F)

final_df <- merge(filtered_df, sbs9_status_df, by = "Sample_name", all.x = TRUE)
final_df <- final_df[!is.na(final_df$SNV_per_callable), ]
final_df <- final_df[!(final_df$Sample_name %in% blacklist), ]

final_df$SBS9_status[is.na(final_df$SBS9_status)] <- "Malignant"

final_df$SBS9_status[final_df$SBS9_status == "Positive"] <- "WT_SBS9_pos"
final_df$SBS9_status[final_df$SBS9_status == "Negative"]  <- "WT_SBS9_neg"

colnames(final_df)[colnames(final_df) == "SNV_per_callable"] <- "load"
colnames(final_df)[colnames(final_df) == "Age_at_sampling_Y"] <- "age"
colnames(final_df)[colnames(final_df) == "SBS9_status"] <- "cell"
final_df <- final_df[, !(colnames(final_df) %in% c("Sample_name", "Myc_translocation_IGV"))]
final_df$Cohort <- "PMC"

# Import Machado et al. data (https://www.nature.com/articles/s41586-022-05072-7)

Machado_table  <- read.table("~/Downloads/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", 
                             header = TRUE, stringsAsFactors = FALSE, sep = "\t")
Machado_B_table <- subset(Machado_table, !(Cell.type2 %in% c("Treg", "HSC", "Naive T", "Memory T")))
Machado_B_table <- subset(Machado_B_table, !(colony %in% c("PD40667rx", "PD40667vu")))

# left with only naive and memory B cells

snvs.list <- round(Machado_B_table$Nmut_hsc_as, digits = 0)
age.list <- Machado_B_table$Age
cell.list <- Machado_B_table$Cell.type2

df.si <- data.frame(age.list, cell.list, snvs.list)
colnames(df.si) <- c("age", "cell", "load")
df.si$Cohort <- "Machado"
df.si$Biopsy_type <- "Healthy_tissue"

# add Bulk PMC

bulk_wgs_meta <- read_xlsx("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_manuscript.xlsx")
bulk_wgs_clonal <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Bulk/Data/clonal_counts_per_sample.csv")
bulk_wgs_meta <- bulk_wgs_meta[, c("PMC_ID", "Age_at_sampling" , "Tumor")]
bulk_merged_df <- bulk_wgs_clonal %>%
  inner_join(bulk_wgs_meta, by = c("Sample" = "Tumor"))

bulk_merged_df <- bulk_merged_df %>%
  filter(!grepl("KOD", Sample))  # outlier

bulk_merged_df <- bulk_merged_df %>%
  dplyr::rename(
    age  = Age_at_sampling,  # or Age_at_sampling_Y if that's the column name
    load = n_clonal,
    cell = PMC_ID
  )
bulk_merged_df$cell <- "Bulk"
bulk_merged_df$Biopsy_type <- "Bulk"
bulk_merged_df$Cohort <- "PMC"
bulk_merged_df$Sample <- NULL

df <- rbind(final_df, df.si, bulk_merged_df)
df$cell_cohort <- paste(df$cell, df$Cohort, sep = "_")

MRCA_load <- read_xlsx("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/MRCA_load.xlsx")

df <- rbind(df, MRCA_load)
df["donor"] <- as.factor(df$age)

# Subset: Naive B from Machado
df.nbc <- subset(df, cell_cohort %in% c("Naive B_Machado"))
df.nbc <- droplevels(df.nbc)

library(nlme)
model_nbc <- lme(fixed = load ~ age, random = ~ 1 | donor, data = df.nbc, method = "ML")

# Extract slope/intercept and (optional) p-value
fx <- fixef(model_nbc)
intercept <- fx["(Intercept)"]
slope     <- fx["age"]
p_val     <- summary(model_nbc)$tTable["age", "p-value"]

# Scatterplot with fitted line
library(ggplot2)
ggplot(df.nbc, aes(x = age, y = load)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_abline(intercept = intercept, slope = slope, linewidth = 1) +
  annotate("text",
           x = min(df.nbc$age, na.rm = TRUE),
           y = max(df.nbc$load, na.rm = TRUE),
           label = paste0("slope = ", round(slope, 1),
                          "\nintercept = ", round(intercept, 1),
                          "\np = ", signif(p_val, 3)),
           hjust = 0, vjust = 1, size = 5) +
  theme_minimal() +
  theme_CHemALL() +
  labs(title = "Naive B cells (Machado)",
       x = "Age (years)", y = "Mutational load") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 16),      # axis tick labels
    axis.title = element_text(size = 18),     # axis titles
    plot.title = element_text(size = 20, face = "bold") # plot title
  )

ggplot(df.nbc, aes(x = age, y = load)) +
  geom_point(size = 3, alpha = 0.85) +
  # add the regression line
  geom_abline(intercept = intercept, slope = slope, linewidth = 1) +
  # add extra points from df where cell_cohort == "WT_SBS9_neg_PMC"
  geom_point(
    data = subset(df, cell_cohort %in% "WT_SBS9_neg_PMC"),
    aes(x = age, y = load),
    color = "red", size = 3, shape = 17
  ) +
  annotate("text",
           x = min(df.nbc$age, na.rm = TRUE),
           y = max(df.nbc$load, na.rm = TRUE),
           label = paste0("slope = ", round(slope, 1),
                          "\nintercept = ", round(intercept, 1),
                          "\np = ", signif(p_val, 3)),
           hjust = 0, vjust = 1, size = 5) +
  theme_minimal() +
  theme_CHemALL() +
  labs(title = "Naive B cells (Machado)",
       x = "Age (years)", y = "Mutational load") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 16),      
    axis.title = element_text(size = 18),     
    plot.title = element_text(size = 20, face = "bold") 
  )

# Subset: Memory B from 'cell_cohort'
df.nbcmem <- subset(df, cell_cohort %in% c("Memory B_Machado"))
df.nbcmem <- droplevels(df.nbcmem)

library(nlme)
model_nbcmem <- lme(fixed = load ~ age, random = ~ 1 | donor, data = df.nbcmem, method = "ML")

# Extract slope/intercept and (optional) p-value
fx_mem <- fixef(model_nbcmem)
intercept_mem <- fx_mem["(Intercept)"]
slope_mem    <- fx_mem["age"]
p_val_mem    <- summary(model_nbcmem)$tTable["age", "p-value"]

# Scatterplot with fitted line
library(ggplot2)
ggplot(df.nbcmem, aes(x = age, y = load)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_abline(intercept = intercept_mem, slope = slope_mem, linewidth = 1) +
  annotate("text",
           x = min(df.nbcmem$age, na.rm = TRUE),
           y = max(df.nbcmem$load, na.rm = TRUE),
           label = paste0("slope = ", round(slope_mem, 1),
                          "\nintercept = ", round(intercept_mem, 1),
                          "\np = ", signif(p_val_mem, 3)),
           hjust = 0, vjust = 1, size = 5) +
  theme_minimal() +
  theme_CHemALL() +
  labs(title = "Memory B cells (Machado)",
       x = "Age (years)", y = "Mutational load") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 16),      # axis tick labels
    axis.title = element_text(size = 18),     # axis titles
    plot.title = element_text(size = 20, face = "bold") # plot title
  )

ggplot(df.nbcmem, aes(x = age, y = load)) +
  geom_point(size = 3, alpha = 0.85) +
  
  # Regression line
  geom_abline(intercept = intercept_mem, slope = slope_mem, linewidth = 1) +
  
  # Add WT_SBS9_pos_PMC in red
  geom_point(
    data = subset(df, cell_cohort == "WT_SBS9_pos_PMC"),
    aes(x = age, y = load),
    color = "red", size = 3, shape = 17
  ) +
  
  # Add MRCA_PMC in blue
  geom_point(
    data = subset(df, cell_cohort == "MRCA_PMC"),
    aes(x = age, y = load),
    color = "blue", size = 3, shape = 17
  ) +
  
  # Add Bulk_PMC in green
  geom_point(
    data = subset(df, cell_cohort == "Bulk_PMC"),
    aes(x = age, y = load),
    color = "green", size = 3, shape = 17
  ) +
  
  annotate("text",
           x = min(df.nbcmem$age, na.rm = TRUE),
           y = max(df.nbcmem$load, na.rm = TRUE),
           label = paste0("slope = ", round(slope_mem, 1),
                          "\nintercept = ", round(intercept_mem, 1),
                          "\np = ", signif(p_val, 3)),
           hjust = 0, vjust = 1, size = 5) +
  
  theme_minimal() +
  theme_CHemALL() +
  labs(title = "Memory B cells (Machado)",
       x = "Age (years)", y = "Mutational load") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold")
  )

ggplot(df.nbcmem, aes(x = age, y = load)) +
  
  # Keep only the regression line (from memory B-cell fit)
  geom_abline(intercept = intercept_mem, slope = slope_mem, linewidth = 1) +
  
  # Add malignant cells as black triangles
  geom_point(
    data = subset(df, cell_cohort == "Malignant_PMC"),
    aes(x = age, y = load),
    color = "black", size = 3, shape = 17
  ) +
  
  # Add WT_SBS9_pos_PMC in red triangles
  geom_point(
    data = subset(df, cell_cohort == "WT_SBS9_pos_PMC"),
    aes(x = age, y = load),
    color = "red", size = 3, shape = 17
  ) +
  
  # Add MRCA_PMC in blue triangles
  geom_point(
    data = subset(df, cell_cohort == "MRCA_PMC"),
    aes(x = age, y = load),
    color = "blue", size = 3, shape = 17
  ) +
  
  # Add Bulk_PMC in green triangles
  geom_point(
    data = subset(df, cell_cohort == "Bulk_PMC"),
    aes(x = age, y = load),
    color = "green", size = 3, shape = 17
  ) +
  
  # Add regression stats
  annotate("text",
           x = min(df.nbcmem$age, na.rm = TRUE),
           y = max(df.nbcmem$load, na.rm = TRUE),
           label = paste0("slope = ", round(slope_mem, 1),
                          "\nintercept = ", round(intercept_mem, 1),
                          "\np = ", signif(p_val, 3)),
           hjust = 0, vjust = 1, size = 5) +
  
  # Set age axis from 0 to 100
  scale_x_continuous(limits = c(0, 100)) +
  
  theme_minimal() +
  theme_CHemALL() +
  labs(title = "Memory B cells (Machado)",
       x = "Age (years)", y = "Mutational load") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold")
  )

# Use the Naive B baseline to compute expected & OE for ALL samples
df$expected_muts_naive <- intercept + slope * df$age
df$OEratio <- df$load / df$expected_muts_naive

df$expected_muts_mem <- intercept_mem + slope_mem * df$age
df$OEratio_mem <- df$load / df$expected_muts_mem

# Build split x variable: Malignant_PMC -> "Malignant_PMC_age<age>"

#df_split <- df %>%
#  filter(!(cell_cohort == "WT_SBS9_neg_PMC" & age == 17.0)) %>%
#  mutate(
#   age_collapsed = case_when(
#      cell_cohort == "Malignant_PMC" & (dplyr::near(age, 14) | dplyr::near(age, 14.7)) ~ 14,
#      TRUE ~ age
#    ),
#    age_lab = format(round(age_collapsed, 1), nsmall = 1),
#    cell_cohort_split = if_else(
#      cell_cohort == "Malignant_PMC",
#      paste0("Malignant_PMC_age", age_lab),
#      cell_cohort
#     )
#  ) %>% droplevels()

df_split <- df %>%
  #filter(!(cell_cohort == "WT_SBS9_neg_PMC" & age == 17.0)) %>%
  mutate(
    age_lab = format(round(age, 1), nsmall = 1),
    cell_cohort_split = case_when(
      cell_cohort == "Malignant_PMC" & dplyr::near(age, 4.1) ~ 
        paste0("Malignant_PMC_age", age_lab, "_", Biopsy_type),
      cell_cohort == "Malignant_PMC" ~ 
        paste0("Malignant_PMC_age", age_lab),
      TRUE ~ cell_cohort
    )
  ) %>%
  droplevels()

# Order factor levels
wt_levels <- c("Naive B_Machado", "WT_SBS9_neg_PMC", "Memory B_Machado", "WT_SBS9_pos_PMC", "Bulk_PMC")

malignant_levels <- df_split %>%
  filter(cell_cohort == "Malignant_PMC") %>%
  arrange(as.numeric(age_lab)) %>%
  pull(cell_cohort_split) %>%
  unique()

df_split <- df_split %>%
  mutate(cell_cohort_split = factor(cell_cohort_split, levels = c(wt_levels, malignant_levels)))

# comparisons
my_comparisons_split <- c(
  list(c("Naive B_Machado", "WT_SBS9_neg_PMC")),
  list(c("Memory B_Machado", "WT_SBS9_pos_PMC")),
  list(c("WT_SBS9_pos_PMC", "Bulk_PMC")),
 # list(c("WT_SBS9_pos_PMC", "MRCA_PMC")),
 # list(c("Bulk_PMC", "MRCA_PMC")),
  lapply(malignant_levels, function(m) c("WT_SBS9_pos_PMC", m))
)

my_colors <- c(
  "Malignant_PMC"    = "#4378bd",
  "Bulk_PMC"    = "grey",
  "MRCA_PMC"    = "grey",
  "WT_SBS9_neg_PMC"  = "#e7872b",
  "WT_SBS9_pos_PMC"  = "#cc6d20", 
  "Naive B_Machado"  = "lightgrey",
  "Memory B_Machado" = "darkgrey"
)

ggplot(
  df_split,
  aes(
    x = factor(cell_cohort_split,
               levels = c(
                 "Naive B_Machado",
                 "WT_SBS9_neg_PMC",
                 "Memory B_Machado",
                 "WT_SBS9_pos_PMC",
                 "Bulk_PMC",
                 "Malignant_PMC_age 4.1_LN",
                 "Malignant_PMC_age 4.1_BM",
                 "Malignant_PMC_age 6.6",
                 "Malignant_PMC_age12.7",
                 "Malignant_PMC_age13.8",
                 "Malignant_PMC_age14.0",
                 "Malignant_PMC_age14.7",
                 "Malignant_PMC_age17.0"
               )),
    y = OEratio,
    #y = OEratio_mem,
    fill = cell_cohort
  )
) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = my_colors) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme_CHemALL() +
  labs(x = "Cell type (Malignant split by age)", y = "OE ratio") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12, angle = 40, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title  = element_text(size = 18, face = "bold", hjust = 0.5)
  ) +
  stat_compare_means(
    comparisons = my_comparisons_split,
    method = "wilcox.test",
    label = "p.format",
    size = 3
  ) +
  scale_y_continuous(breaks = seq(0, 25, 1), limits = c(0, 25))

median_table <- df_split %>%
  group_by(cell_cohort_split) %>%
  summarise(median_OE = median(OEratio, na.rm = TRUE))


##########

# Keep only malignant cells and clean
df_mal <- df %>%
  filter(cell_cohort == "Malignant_PMC") %>%
  filter(!is.na(OEratio), !is.na(age)) %>%
  droplevels()


ggplot(df_mal, aes(x = age, y = load)) +
  geom_point(size = 3, alpha = 0.85) +
  theme_minimal() +
  theme_CHemALL() +
  labs(title = "Malignant B cells (PMC)",
       x = "Age (years)", y = "Mutational load") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 16),      # axis tick labels
    axis.title = element_text(size = 18),     # axis titles
    plot.title = element_text(size = 20, face = "bold") # plot title
  )

# Values to subtract per age
subtract_values <- tibble(
  age = c(4.1, 6.6, 12.7, 13.8, 14.0,14.7, 17.0),
  subtract_load = c(687, 725, 487, 1206, 696,949, 354)
)

# Filter malignant cells and subtract age-specific MRCA load
df_mal_without_MRCA <- df %>%
  filter(cell_cohort == "Malignant_PMC") %>%
  filter(!is.na(load), !is.na(age)) %>%
  left_join(subtract_values, by = "age") %>%
  mutate(
    subtract_load = ifelse(is.na(subtract_load), 0, subtract_load),  # if age not in list → subtract 0
    load_corrected = load - subtract_load
  ) %>%
  droplevels()

# Plot corrected mutational load with jitter
ggplot(df_mal_without_MRCA, aes(x = age, y = load_corrected)) +
  geom_jitter(width = 0.2, height = 0, size = 3, alpha = 0.85) +   # jitter added here
  theme_minimal() +
  theme_CHemALL() +
  labs(title = "Malignant B cells (PMC)",
       x = "Age (years)", 
       y = "Mutational load (corrected)") +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold")
  )


# LME: O/E ~ age with random intercept per donor
model_mal <- nlme::lme(
  OEratio ~ age,
  random = ~ 1 | donor,          # use Patient if that's your grouping column
  data   = df_mal,
  method = "REML",
  na.action = na.omit
)

# Fixed effects + p-value
fx     <- fixef(model_mal)
b0     <- unname(fx["(Intercept)"])
b1     <- unname(fx["age"])
p_age  <- summary(model_mal)$tTable["age", "p-value"]

# Population-level predictions + 95% CI from fixed-effects vcov
newdat <- data.frame(age = seq(min(df_mal$age), max(df_mal$age), length.out = 100))
# Fixed-effect predictions + SE from vcov
fx <- fixef(model_mal)
V  <- vcov(model_mal)
X  <- model.matrix(~ age, data = newdat)

newdat$fit <- as.numeric(X %*% fx)
se <- sqrt(diag(X %*% V %*% t(X)))
newdat$lwr <- newdat$fit - 1.96 * se
newdat$upr <- newdat$fit + 1.96 * se
# Scatter + LME fit
ggplot(df_mal, aes(x = age, y = OEratio)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2.8, alpha = 0.9) +
  geom_ribbon(data = newdat,
              aes(x = age, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.15) +
  geom_line(data = newdat, aes(y = fit), linewidth = 1) +
  theme_minimal() + theme_CHemALL() +
  labs(title = "Malignant PMC — O/E ratio vs age (LME)",
       x = "Age (years)", y = "O/E ratio") +
  annotate("text",
           x = 15,
           y = max(df_mal$OEratio, na.rm = TRUE) * 0.98,
           label = paste0("slope = ", round(b1, 3),
                          "\nintercept = ", round(b0, 3),
                          "\np = ", signif(p_age, 3)),
           hjust = 0, vjust = 1, size = 5) +
  expand_limits(x = 15) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  theme(
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    plot.title   = element_text(size = 20, face = "bold")
  )


# # Subset: SBS9 neg WT B from PMC
# df.nbcsbs9neg <- subset(df, cell_cohort %in% c("WT_SBS9_neg_PMC"))
# df.nbcsbs9neg <- droplevels(df.nbcsbs9neg)
# 
# model_nbcsbs9neg <- lme(fixed = load ~ age, random = ~ 1 | donor, data = df.nbcsbs9neg, method = "ML")
# 
# # Extract slope/intercept and (optional) p-value
# fx_sbs9neg <- fixef(model_nbcsbs9neg)
# intercept_sbs9neg <- fx_sbs9neg["(Intercept)"]
# slope_sbs9neg   <- fx_sbs9neg["age"]
# p_val_sbs9neg   <- summary(model_nbcsbs9neg)$tTable["age", "p-value"]
# 
# # Scatterplot with fitted line
# ggplot(df.nbcsbs9neg, aes(x = age, y = load)) +
#   geom_point(size = 3, alpha = 0.85) +
#   geom_abline(intercept = intercept_sbs9neg, slope = slope_sbs9neg, linewidth = 1) +
#   annotate("text",
#            x = min(df.nbcsbs9neg$age, na.rm = TRUE),
#            y = max(df.nbcsbs9neg$load, na.rm = TRUE),
#            label = paste0("slope = ", round(slope_sbs9neg, 1),
#                           "\nintercept = ", round(intercept_sbs9neg, 1),
#                           "\np = ", signif(p_val_sbs9neg, 3)),
#            hjust = 0, vjust = 1, size = 5) +
#   theme_minimal() +
#   theme_CHemALL() +
#   labs(title = "SBS9neg WT B cells (Burkitt)",
#        x = "Age (years)", y = "Mutational load") +
#   theme(
#     legend.position = "none",
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.text = element_text(size = 16),      # axis tick labels
#     axis.title = element_text(size = 18),     # axis titles
#     plot.title = element_text(size = 20, face = "bold") # plot title
#   )
# 
# # Subset: SBS9 neg WT B from PMC without PIA9 
# df.nbcsbs9neg_wo_PIA9 <- subset(df, cell_cohort %in% c("WT_SBS9_neg_PMC") & age != 17)
# df.nbcsbs9neg_wo_PIA9 <- droplevels(df.nbcsbs9neg_wo_PIA9)
# 
# model_nbcsbs9neg_wo_PIA9 <- lme(fixed = load ~ age, random = ~ 1 | donor, data = df.nbcsbs9neg_wo_PIA9, method = "ML")
# 
# # Extract slope/intercept and (optional) p-value
# fx_sbs9neg_wo_PIA9 <- fixef(model_nbcsbs9neg_wo_PIA9)
# intercept_sbs9neg_wo_PIA9 <- fx_sbs9neg_wo_PIA9["(Intercept)"]
# slope_sbs9neg_wo_PIA9   <- fx_sbs9neg_wo_PIA9["age"]
# p_val_sbs9neg_wo_PIA9   <- summary(model_nbcsbs9neg_wo_PIA9)$tTable["age", "p-value"]
# 
# # Scatterplot with fitted line
# ggplot(df.nbcsbs9neg_wo_PIA9, aes(x = age, y = load)) +
#   geom_point(size = 3, alpha = 0.85) +
#   geom_abline(intercept = intercept_sbs9neg_wo_PIA9, slope = slope_sbs9neg_wo_PIA9, linewidth = 1) +
#   annotate("text",
#            x = min(df.nbcsbs9neg_wo_PIA9$age, na.rm = TRUE),
#            y = max(df.nbcsbs9neg_wo_PIA9$load, na.rm = TRUE),
#            label = paste0("slope = ", round(slope_sbs9neg_wo_PIA9, 1),
#                           "\nintercept = ", round(intercept_sbs9neg_wo_PIA9, 1),
#                           "\np = ", signif(p_val_sbs9neg, 3)),
#            hjust = 0, vjust = 1, size = 5) +
#   theme_minimal() +
#   theme_CHemALL() +
#   labs(title = "SBS9neg WT B cells (Burkitt)_wo_PIA9",
#        x = "Age (years)", y = "Mutational load") +
#   theme(
#     legend.position = "none",
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.text = element_text(size = 16),      # axis tick labels
#     axis.title = element_text(size = 18),     # axis titles
#     plot.title = element_text(size = 20, face = "bold") # plot title
#   )
# 
# # Subset: SBS9 pos WT B from PMC
# df.nbcsbs9pos <- subset(df, cell_cohort %in% c("WT_SBS9_pos_PMC"))
# df.nbcsbs9pos <- droplevels(df.nbcsbs9pos)
# 
# model_nbcsbs9pos <- lme(fixed = load ~ age, random = ~ 1 | donor, data = df.nbcsbs9pos, method = "ML")
# 
# # Extract slope/intercept and (optional) p-value
# fx_sbs9pos <- fixef(model_nbcsbs9pos)
# intercept_sbs9pos <- fx_sbs9pos["(Intercept)"]
# slope_sbs9pos   <- fx_sbs9pos["age"]
# p_val_sbs9pos   <- summary(model_nbcsbs9pos)$tTable["age", "p-value"]
# 
# # Scatterplot with fitted line
# ggplot(df.nbcsbs9pos, aes(x = age, y = load)) +
#   geom_point(size = 3, alpha = 0.85) +
#   geom_abline(intercept = intercept_sbs9pos, slope = slope_sbs9pos, linewidth = 1) +
#   annotate("text",
#            x = min(df.nbcsbs9pos$age, na.rm = TRUE),
#            y = max(df.nbcsbs9pos$load, na.rm = TRUE),
#            label = paste0("slope = ", round(slope_sbs9pos, 1),
#                           "\nintercept = ", round(intercept_sbs9pos, 1),
#                           "\np = ", signif(p_val_sbs9pos, 3)),
#            hjust = 0, vjust = 1, size = 5) +
#   theme_minimal() +
#   theme_CHemALL() +
#   labs(title = "SBS9pos WT B cells (Burkitt)",
#        x = "Age (years)", y = "Mutational load") +
#   theme(
#     legend.position = "none",
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.text = element_text(size = 16),      # axis tick labels
#     axis.title = element_text(size = 18),     # axis titles
#     plot.title = element_text(size = 20, face = "bold") # plot title
#   )
# 
# # Subset: malignant from PMC
# df.mal <- subset(df, cell_cohort %in% c("Malignant_PMC"))
# df.mal <- droplevels(df.mal)
# 
# model_mal <- lme(fixed = load ~ age, random = ~ 1 | donor, data = df.mal, method = "ML")
# 
# # Extract slope/intercept and (optional) p-value
# fx_mal <- fixef(model_mal)
# intercept_mal <- fx_mal["(Intercept)"]
# slope_mal   <- fx_mal["age"]
# p_val_mal   <- summary(model_mal)$tTable["age", "p-value"]
# 
# # Scatterplot with fitted line
# ggplot(df.mal, aes(x = age, y = load)) +
#   geom_point(size = 3, alpha = 0.85) +
#   geom_abline(intercept = intercept_mal, slope = slope_mal, linewidth = 1) +
#   annotate("text",
#            x = min(df.mal$age, na.rm = TRUE),
#            y = max(df.mal$load, na.rm = TRUE),
#            label = paste0("slope = ", round(slope_mal, 1),
#                           "\nintercept = ", round(intercept_mal, 1),
#                           "\np = ", signif(p_val_mal, 3)),
#            hjust = 0, vjust = 1, size = 5) +
#   theme_minimal() +
#   theme_CHemALL() +
#   labs(title = "Malignant B cells (Burkitt)",
#        x = "Age (years)", y = "Mutational load") +
#   theme(
#     legend.position = "none",
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.text = element_text(size = 16),      # axis tick labels
#     axis.title = element_text(size = 18),     # axis titles
#     plot.title = element_text(size = 20, face = "bold") # plot title
#   )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###
# 
# # Use this SBS9neg WT B cell from Burkitts as a baseline to compute expected & OE for all Burkitt samples
# df$expected_muts_SBS9neg <- intercept_sbs9neg + slope_sbs9neg * df$age
# df$OEratio_SBS9neg <- df$load / df$expected_muts_SBS9neg
# 
# # comparisons must exactly match levels(df$cell)
# my_comparisons_new <- list(
#   c("WT_SBS9_neg_PMC", "WT_SBS9_pos_PMC"),
#   c("WT_SBS9_pos_PMC", "Malignant_PMC"))
# 
# ggplot(
#   subset(df, grepl("PMC", Cohort)),
#   aes(x = cell_cohort, y = OEratio_SBS9neg, fill = cell_cohort)
# ) +
#   geom_boxplot() +
#   scale_fill_manual(values = my_colors) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   theme_CHemALL() +
#   labs(x = "Cell type", y = "OE ratio") +
#   theme(
#     legend.position = "none",
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.text.x = element_text(size = 14),
#     axis.text.y = element_text(size = 14),
#     axis.title.x = element_text(size = 16, face = "bold"),
#     axis.title.y = element_text(size = 16, face = "bold"),
#     plot.title  = element_text(size = 18, face = "bold", hjust = 0.5)
#   ) +
#   stat_compare_means(
#     comparisons = my_comparisons_new,
#     method = "wilcox.test",
#     label = "p.format",
#     size = 5  ) +
#   scale_y_continuous(breaks = seq(0, 15, 1), limits = c(0, 15))
# 
# # Keep PMC rows and make a split x-variable only for Malignant_PMC
# df_pmc <- df %>%
#   filter(grepl("PMC", Cohort)) %>%
#   mutate(
#     # collapse 14 and 14.7 into the same age label "14"
#     age_collapsed = case_when(
#       cell_cohort == "Malignant_PMC" & (dplyr::near(age, 14) | dplyr::near(age, 14.7)) ~ 14,
#       TRUE ~ age
#     ),
#     # label includes the collapsed age (no donor in label per your current code)
#     cell_cohort_donor = if_else(
#       cell_cohort == "Malignant_PMC",
#       paste0("Malignant_PMC_", "age", age_collapsed),
#       cell_cohort
#     )
#   ) %>%
#   filter(!is.na(OEratio_SBS9neg)) %>%
#   droplevels()
# 
# # Order by (collapsed) age for Malignant_PMC, keep WT_* first
# wt_levels <- c("WT_SBS9_neg_PMC", "WT_SBS9_pos_PMC")
# 
# malignant_levels <- df_pmc %>%
#   filter(cell_cohort == "Malignant_PMC") %>%
#   arrange(age_collapsed) %>%
#   pull(cell_cohort_donor) %>%
#   unique()
# 
# df_pmc <- df_pmc %>%
#   mutate(cell_cohort_donor = factor(cell_cohort_donor,
#                                     levels = c(wt_levels, malignant_levels)))
# 
# # Add WT vs WT, then all WT_pos vs each malignant-age box
# my_comparisons_split <- c(
#   list(c("WT_SBS9_neg_PMC", "WT_SBS9_pos_PMC")),   # <-- WT comparison
#   lapply(malignant_levels, function(m) c("WT_SBS9_pos_PMC", m))
# )
# 
# # Plot
# ggplot(
#   df_pmc,
#   aes(x = cell_cohort_donor, y = OEratio_SBS9neg, fill = cell_cohort)
# ) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +
#   geom_point(position = position_jitter(width = 0.15, height = 0),
#              size = 2.8, alpha = 0.9, stroke = 0.4) +
#   scale_fill_manual(values = my_colors) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   theme_CHemALL() +
#   labs(x = "Cell type / donor (Malignant split by age)", y = "OE ratio") +
#   theme(
#     legend.position = "right",
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.text.x = element_text(size = 12, angle = 40, hjust = 1),
#     axis.text.y = element_text(size = 14),
#     axis.title.x = element_text(size = 16, face = "bold"),
#     axis.title.y = element_text(size = 16, face = "bold"),
#     plot.title  = element_text(size = 18, face = "bold", hjust = 0.5)
#   ) +
#   stat_compare_means(
#     comparisons = my_comparisons_split,
#     method = "wilcox.test",
#     label = "p.format",
#     size = 5
#   ) +
#   coord_cartesian(ylim = c(0, 15))
# 
# ggplot(
#   subset(df, grepl("PMC", Cohort) & !(cell_cohort == "WT_SBS9_neg_PMC" & age == 17)),
#   aes(x = cell_cohort, y = OEratio_SBS9neg, fill = cell_cohort)
#   ) +
#   geom_boxplot() +
#   scale_fill_manual(values = my_colors) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   theme_CHemALL() +
#   labs(x = "Cell type", y = "OE ratio") +
#   theme(
#     legend.position = "none",
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.text.x = element_text(size = 14),
#     axis.text.y = element_text(size = 14),
#     axis.title.x = element_text(size = 16, face = "bold"),
#     axis.title.y = element_text(size = 16, face = "bold"),
#     plot.title  = element_text(size = 18, face = "bold", hjust = 0.5)
#   ) +
#   stat_compare_means(
#     comparisons = my_comparisons_new,
#     method = "wilcox.test",
#     label = "p.format",
#     size = 5  ) +
#   scale_y_continuous(breaks = seq(0, 15, 1), limits = c(0, 15))
# 
# 
# ###
# # Bulk
# 
# bulk_wgs_meta <- read_xlsx("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_manuscript.xlsx")
# bulk_wgs_clonal <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Bulk/Data/clonal_counts_per_sample.csv")
# bulk_wgs_meta <- bulk_wgs_meta[, c("PMC_ID", "Age_at_sampling" , "Tumor")]
# merged_df <- bulk_wgs_clonal %>%
#   inner_join(bulk_wgs_meta, by = c("Sample" = "Tumor"))
# 
# merged_df <- merged_df %>%
#   filter(!grepl("KOD", Sample))  # outlier
# 
# merged_df <- merged_df %>%
#   dplyr::rename(
#     age  = Age_at_sampling,  # or Age_at_sampling_Y if that's the column name
#     load = n_clonal,
#     cell = PMC_ID
#   )
# merged_df$cell <- "Bulk"
# 
# merged_df$Cohort <- "PMC"
# 
# df_all <- bind_rows(df, merged_df)
# df_all$cell_cohort <- paste(df_all$cell, df_all$Cohort, sep = "_")
# 
# df_all["donor"] <- as.factor(df_all$age)
# 
# df_all$expected_muts_naive <- intercept + slope * df_all$age
# df_all$OEratio <- df_all$load / df_all$expected_muts_naive
# 
# 
# 
# ggplot(df_all, aes(x = cell_cohort, y = OEratio, fill = cell_cohort)) +
#   geom_boxplot() +
#   scale_fill_manual(values = my_colors) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
#   theme_minimal() +
#   theme_CHemALL() +
#   labs(x = "Cell type", y = "OE ratio") +
#   theme(
#     legend.position = "none",
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank(),
#     axis.text.x = element_text(size = 14),
#     axis.text.y = element_text(size = 14),
#     axis.title.x = element_text(size = 16, face = "bold"),
#     axis.title.y = element_text(size = 16, face = "bold"),
#     plot.title  = element_text(size = 18, face = "bold", hjust = 0.5)
#   ) +
#   #stat_compare_means(
#   #  comparisons = my_comparisons,
#   #  method = "wilcox.test",
#   #  label = "p.format",
#   #  size = 5  ) +
#   scale_y_continuous(breaks = seq(0, 14, 1), limits = c(0, 14))
# 
# 
