################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to plot oncoplot for P3G6
# Author: Alexander Steemers
# Date: June 2025
################################################################################

# Load libraries

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(readxl)

# Set working directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/")

# Load drivers found using the script: /hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Driver_analysis_P3G6.sh

drivers_P3G6 <- read.csv("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/Drivers/P3G6_drivers.csv")
drivers_P3G6 <- unique(drivers_P3G6)

# Export full driver list

write.csv(drivers_P3G6, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/P3G6_fulldriverlist.csv", row.names = F)

# Load all known Burkitt lymphoma driver genes from IntOGen #https://www.intogen.org/search?cancer=BL

intogen_driver_list <- list.files(
  path = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Drivers",
  pattern = "\\.tsv$",
  full.names = TRUE
)

intogen_driver_list <- grep("DriverGenes_BL\\.tsv$", intogen_driver_list, value = TRUE)

# Read and combine all TSV files into one data frame

drivers_intogen <- do.call(rbind, lapply(intogen_driver_list, function(file) {
  read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}))

# Load known BL driver genes from The whole-genome landscape of Burkitt lymphoma subtypes #https://ashpublications.org/blood/article/134/19/1598/375002/The-whole-genome-landscape-of-Burkitt-lymphoma

blood2019_drivers <- read_xlsx("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Drivers/bloodbld2019001880-suppl2.xlsx", sheet = 5, skip = 2)

# Merge these drivers from two datasets

drivers_burkitt <- unique(c(drivers_intogen$Symbol, blood2019_drivers$Genes))

# Filter drivers based on driver list

drivers_P3G6_f <- drivers_P3G6[drivers_P3G6$Gene %in% drivers_burkitt, ]

drivers_P3G6_dedup <- drivers_P3G6_f[!duplicated(drivers_P3G6_f[c("Sample", "Gene", "Chrom", "Pos", "Ref", "Alt")]), ]

# Change names of effects

drivers_P3G6_shared <- drivers_P3G6_dedup %>%
  mutate(Effect = case_when(
    Effect == "missense_variant"    ~ "missense",
    Effect == "missense_variant&splice_region_variant"    ~ "missense",
    Effect == "stop_gained"         ~ "nonsense",
    Effect == "inframe_insertion"   ~ "indel",
    Effect == "frameshift_variant"   ~ "indel",
    Effect == "downstream_gene_variant"    ~ "missense",
    TRUE                            ~ Effect  # keep unchanged if not matched
  ))

# Create the matrix with effects as values, and combine effects with ';' if needed

effect_matrix <- drivers_P3G6_shared %>%
  dplyr::select(Gene, Sample, Effect) %>%
  distinct() %>%
  group_by(Gene, Sample) %>%
  summarise(Effect = paste(Effect, collapse = "; "), .groups = "drop") %>%
  pivot_wider(
    names_from = Gene,  # Gene names as columns
    values_from = Effect,
    values_fill = list(Effect = " ")  # Fill missing values with empty string
  )

# Convert tibble to data frame

effect_matrix <- as.data.frame(effect_matrix)

# Set Sample as row names

rownames(effect_matrix) <- effect_matrix$Sample

# Drop the Sample column as it's now the row names

effect_matrix <- effect_matrix[, -1]
effect_matrix_t <- t(as.matrix(effect_matrix))

# Give colours and band heights to the different mutation types

col = c("indel" = "blue", "nonsense" = "red", "missense" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # small blue
  indel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["indel"], col = NA))
  },
  # big red
  nonsense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["nonsense"], col = NA))
  },
  # small green
  missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["missense"], col = NA))
  }
)

# Provide sample info

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
P3G6_rows <- input_df %>% filter(Novogene_ID == "P3G6") %>% arrange(factor(Myc_translocation_IGV, levels = c("no", "yes")))   
P3G6_rows <- P3G6_rows %>% filter(!grepl("MSCBULK", Sample_name))
P3G6_samples <- P3G6_rows$Sample_name

# Filter the sample vector to include only those present in the matrix

ordered_samples <- P3G6_samples[P3G6_samples %in% colnames(effect_matrix_t)]

# Reorder the columns of the matrix

effect_matrix_t <- effect_matrix_t[, ordered_samples]

# Define annotation colors

group_colors <- c(Ascites = "black")  
sample_group <- ifelse(ordered_samples %in% P3G6_samples, "Ascites", NA)

# Create the top annotation with visible legend

top_annot <- HeatmapAnnotation(
  Location = sample_group,
  col = list(Location = group_colors),
  annotation_name_side = "left"
)

column_title = "P3G6 oncoplot"
heatmap_legend_param = list(title = "Alternations", at = c("indel", "nonsense", "missense"), 
                            labels = c("Indel", "Nonsense", "Missense"))

# Make initial oncoplot

pdf("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/P3G6_oncoprint_before_manual_check.pdf", width = 10, height = 6)  # adjust size as needed
oncoPrint(effect_matrix_t,
          alter_fun = alter_fun,
          col = col,
          column_title = column_title,
          heatmap_legend_param = heatmap_legend_param,
          show_column_names = TRUE,
          column_order = colnames(effect_matrix_t),
          top_annotation = top_annot)
dev.off()

# Make a table of driver mutations to check manually on IGV

drivers_to_check <- drivers_P3G6_shared %>%
  dplyr::select(-Sample) %>%
  distinct()

# Import csv file AFTER manual checking on IGV

effect_matrix_after_check <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/P3G6_manual_driver_check.csv")
effect_matrix_after_check_tibble <- as_tibble(effect_matrix_after_check)

effect_matrix_after_check <- as.data.frame(effect_matrix_after_check_tibble)

# Set Sample as row names

rownames(effect_matrix_after_check) <- effect_matrix_after_check$Sample

# Drop the Sample column as it's now the row names

effect_matrix_after_check <- effect_matrix_after_check[, -1]
effect_matrix_after_check_T <- t(as.matrix(effect_matrix_after_check))

# Convert to character if necessary

colnames(effect_matrix_after_check_T) <- as.character(colnames(effect_matrix_after_check_T))

# Find which sample IDs are missing from the matrix

missing_samples <- setdiff(P3G6_samples, colnames(effect_matrix_after_check_T))

# If there are missing samples, add them as new columns with "" values

new_cols <- matrix("", nrow = nrow(effect_matrix_after_check_T), ncol = length(missing_samples))
colnames(new_cols) <- missing_samples
rownames(new_cols) <- rownames(effect_matrix_after_check_T)
effect_matrix_after_check_T <- cbind(effect_matrix_after_check_T, new_cols)

# Filter the sample vector to include only those present in the matrix

ordered_samples <- P3G6_samples[P3G6_samples %in% colnames(effect_matrix_after_check_T)]

# Define group colors

group_colors <- c(Ascites = "#000000")  

# Define cell type colors

celltype_colors <- c("Normal B-cell" = "#E7872B", "Burkitt Lymphoma cell" = "#3F78C1", "Bulk Tumour" = "#00008B")

# Define sample groups

sample_group <- ifelse(ordered_samples %in% P3G6_samples, "Ascites", NA)

# Define which samples are Normal B cell, which are BL and which are Bulk

normal_bcell_samples <- c("P3G6GPDABC31", "PB11197-BLASC-BCELLP2F4", 
                          "PB11197-BLASC-BCELLP2B4", "PB11197-BLASC-BCELLP2E4", 
                          "PB11197-BLASC-BCELLP2D4", "PB11197-BLASC-BCELLP2C4")


P3G6_rows <- P3G6_rows %>% 
  mutate(
    Category = case_when(
      Sample_name %in% normal_bcell_samples      ~ "Normal B-cell",
      grepl("BULK", Sample_name, ignore.case = TRUE) ~ "Bulk Tumour", 
      TRUE                                        ~ "Burkitt Lymphoma cell"
    )
  )

P3G6_rows <- P3G6_rows %>% 
  mutate(
    Category = factor(Category,
                      levels = c("Normal B-cell",
                                 "Burkitt Lymphoma cell",
                                 "Bulk Tumour"))
  ) %>% 
  arrange(Category)

# Define CellType annotation

cell_type <- as.character(P3G6_rows$Category)

# Create the top annotation with both Location and CellType

top_annot <- HeatmapAnnotation(
  Location = sample_group,
  Cell_Type = cell_type,
  col = list(
    Location = group_colors,
    Cell_Type = celltype_colors
  ),
  annotation_name_side = "left"
)

column_title = "P3G6 oncoplot"
heatmap_legend_param = list(title = "Alternations", at = c("indel", "nonsense", "missense"), 
                            labels = c("Indel", "Nonsense", "Missense"))

ordered_samples <- P3G6_rows$Sample_name

# Re-order the columns of the matrix

effect_matrix_after_check_T <- effect_matrix_after_check_T[ , ordered_samples]

# Print final oncoplot 

pdf("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/P3G6_oncoprint_after_manual_check.pdf", width = 10, height = 6)  # adjust size as needed
oncoPrint(effect_matrix_after_check_T,
          alter_fun = alter_fun,
          col = col,
          column_title = column_title,
          heatmap_legend_param = heatmap_legend_param,
          show_column_names = TRUE,
          column_order = colnames(effect_matrix_after_check_T),
          top_annotation = top_annot,
          right_annotation = NULL)
dev.off()