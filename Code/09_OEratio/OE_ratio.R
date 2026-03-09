################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Script to do O/E ratio
# Author: Alexander Steemers
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
source('Data/theme_burkitt.R')

# Load metadata 
input_df <-  read_excel('Data/Sample_overview.xlsx') 
diagnostic_df <- read.csv('Data/Bulk_sample_manuscript.csv')
low_callable_df     <- "low_callable_loci.csv" # generated after running 03_QC_step1.R
below_curve_df      <- "below_curve_samples.csv" # generated after running 03_QC_step1.R
bad_baf_df          <- "bad_baf_samples.csv" # generated after running 03_QC_step1.R
fail_vaf_df         <- "PTA_samples_failVAFcheck.txt" # generated after running 04_QC_step2.R

# Make blacklist

blacklist <- unique(c(below_curve_df$Sample_name,
                      low_callable_df$Sample_name,
                      bad_baf_df$Sample_name,
                      fail_vaf_df$samplename))

# Load single cell VCF files (PTATO filtered)
# For all_filtered_vcfs take all files listed in Data/single_cell_ptato_filtered_file_names_SNVs.txt which can be in turn found in Mendeley doi: 10.17632/phk5jhhm7d.1 

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

input_df <-  read_excel('Data/Sample_overview.xlsx') 

# Merge with input_df by Sample
merged_df <- merge(input_df, snv_df, by = "Sample_name", all.x = TRUE)

# Normalize counts by callable fraction
merged_df$SNV_per_callable <- merged_df$SNV_count / merged_df$Callable_fraction

filtered_df <- merged_df[, c("Sample_name", "Age_at_sampling_Y", "SNV_per_callable", "Myc_translocation_IGV", "Biopsy_type" )]

sbs9_status_df <- read.csv("Normal_SBS9_status_table.csv", stringsAsFactors =  F) # generated after running wt_pre_post_comparison.R

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

Machado_table  <- read.table("`Data/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001.txt", 
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

bulk_wgs_meta <- read_xlsx("Data/Bulk_sample_manuscript.xlsx")
bulk_wgs_clonal <- read.csv("clonal_counts_per_sample_ccube.csv") # generated after running 05_CFF_bulk.R script

bulk_wgs_meta <- bulk_wgs_meta[, c("Tumor_Sample_Barcode", "Age_at_sampling" , "Tumor")]

bulk_merged_df <- bulk_wgs_clonal %>%
  inner_join(bulk_wgs_meta, by = c("Sample" = "Tumor_Sample_Barcode"))

 bulk_merged_df <- bulk_merged_df %>%
   filter(!grepl("PC1A", Sample))  # outlier

bulk_merged_df <- bulk_merged_df %>%
  dplyr::rename(
    age  = Age_at_sampling,  # or Age_at_sampling_Y if that's the column name
    load = n_clonal,
    cell = Sample
  )
bulk_merged_df$cell <- "Bulk"
bulk_merged_df$Biopsy_type <- "Bulk"
bulk_merged_df$Cohort <- "PMC"
bulk_merged_df$Sample <- NULL
bulk_merged_df <- bulk_merged_df %>%
  select(-Tumor)

df <- rbind(final_df, df.si, bulk_merged_df)
df$cell_cohort <- paste(df$cell, df$Cohort, sep = "_")

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
  theme_burkitt() +
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
  theme_burkitt() +
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
  theme_burkitt() +
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
  theme_burkitt() +
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
  theme_burkitt() +
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
  theme_burkitt() +
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
    size = 3,
    label.y = c(4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)  # adjust numbers
  ) +
  scale_y_continuous(breaks = seq(0, 30, 1), limits = c(0, 30))

median_table <- df_split %>%
  group_by(cell_cohort_split) %>%
  summarise(median_OE = median(OEratio, na.rm = TRUE))


# Now want to split the PIA9 SBS9- normal cells

# comparisons
my_comparisons_split <- c(
  list(c("Naive B_Machado", "WT_SBS9_neg_PMC_with17")),
  list(c("Naive B_Machado", "WT_SBS9_neg_PMC_no17")),
  list(c("WT_SBS9_neg_PMC_no17", "WT_SBS9_pos_PMC")),
  list(c("WT_SBS9_neg_PMC_with17", "WT_SBS9_pos_PMC")),
  list(c("Memory B_Machado", "WT_SBS9_pos_PMC")),
  list(c("WT_SBS9_pos_PMC", "Bulk_PMC")),
  # list(c("WT_SBS9_pos_PMC", "MRCA_PMC")),
  # list(c("Bulk_PMC", "MRCA_PMC")),
  lapply(malignant_levels, function(m) c("WT_SBS9_pos_PMC", m))
)
  
df_plot <- df_split %>%
  mutate(cell_cohort_split2 = cell_cohort_split) %>%
  filter(cell_cohort_split != "WT_SBS9_neg_PMC")

# take only WT_SBS9_neg_PMC rows
wt <- df_split %>%
  filter(cell_cohort_split == "WT_SBS9_neg_PMC")

# version WITH 17.0 (all rows) (PIA9 has age = 17.0)
wt_with17 <- wt %>%
  mutate(cell_cohort_split2 = "WT_SBS9_neg_PMC_with17")

# version WITHOUT 17.0 (drop age 17.0)
wt_no17 <- wt %>%
  filter(age != 17.0) %>%
  mutate(cell_cohort_split2 = "WT_SBS9_neg_PMC_no17")

# Combine everything
df_plot <- bind_rows(
  df_plot,
  wt_with17,
  wt_no17
)

ggplot(
  df_plot,
  aes(
    x = factor(
      cell_cohort_split2,
      levels = c(
        "Naive B_Machado",
        "WT_SBS9_neg_PMC_with17",
        "WT_SBS9_neg_PMC_no17",
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
      )
    ),
    y = OEratio,
    fill = cell_cohort   # same colors as before
  )
) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = my_colors) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  theme_minimal() +
  theme_burkitt() +
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

median_table <- df_plot %>%
  group_by(cell_cohort_split2) %>%
  summarise(median_OE = median(OEratio, na.rm = TRUE))

iqr_table <- df_plot %>%
  group_by(cell_cohort_split2) %>%
  summarise(
    Q1  = quantile(OEratio, 0.25, na.rm = TRUE),
    Q3  = quantile(OEratio, 0.75, na.rm = TRUE),
    IQR = IQR(OEratio, na.rm = TRUE)
  )

##########

# Keep only malignant cells and clean
df_mal <- df %>%
  filter(cell_cohort == "Malignant_PMC") %>%
  filter(!is.na(OEratio), !is.na(age)) %>%
  droplevels()


ggplot(df_mal, aes(x = age, y = load)) +
  geom_point(size = 3, alpha = 0.85) +
  theme_minimal() +
  theme_burkitt() +
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
  theme_burkitt() +
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
  theme_minimal() + theme_burkitt() +
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
