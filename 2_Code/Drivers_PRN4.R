################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to plot oncoplot for PRN4
# Author: Alexander Steemers
# Date: July 2025
################################################################################

# Load libraries

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(readxl)

# Set working directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/")

# Load drivers found using the script: /hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Driver_analysis_PRN4.sh

drivers_PRN4 <- read.csv("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/Drivers/PRN4_drivers.csv")
drivers_PRN4 <- unique(drivers_PRN4)

# Export full driver list

write.csv(drivers_PRN4, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/PRN4_fulldriverlist.csv", row.names = F)

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

drivers_PRN4_f <- drivers_PRN4[drivers_PRN4$Gene %in% drivers_burkitt, ]

drivers_PRN4_dedup <- drivers_PRN4_f[!duplicated(drivers_PRN4_f[c("Sample", "Gene", "Chrom", "Pos", "Ref", "Alt")]), ]

# Change names of effects

drivers_PRN4_shared <- drivers_PRN4_dedup %>%
  mutate(Effect = case_when(
    Effect == "missense_variant"    ~ "missense",
    Effect == "5_prime_UTR_variant"    ~ "missense",
    Effect == "3_prime_UTR_variant&NMD_transcript_variant"    ~ "indel",
    Effect == "stop_gained"         ~ "nonsense",
    Effect == "splice_donor_variant"   ~ "missense",
    Effect == "frameshift_variant"   ~ "indel",
    Effect == "downstream_gene_variant"    ~ "indel",
    Effect == "inframe_deletion&splice_region_variant"    ~ "indel",
    TRUE                            ~ Effect  # keep unchanged if not matched
  ))

# Create the matrix with effects as values, and combine effects with ';' if needed

effect_matrix <- drivers_PRN4_shared %>%
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
PRN4_rows <- input_df %>% filter(Novogene_ID == "PRN4") %>% arrange(factor(Myc_translocation_IGV, levels = c("no", "yes")))   
PRN4_rows <- PRN4_rows %>% filter(!grepl("MSCBULK", Sample_name))
PRN4_samples <- PRN4_rows$Sample_name

PRN4_BM <- PRN4_rows$Sample_name[PRN4_rows$Biopsy_type == "BM"]

PRN4_LN <- PRN4_rows$Sample_name[PRN4_rows$Biopsy_type == "LN"]

# Assign sample group labels

ordered_samples <- PRN4_samples[PRN4_samples %in% colnames(effect_matrix_t)]
all_samples <- colnames(effect_matrix_t)
sample_group <- ifelse(ordered_samples %in% PRN4_LN, "LN",
                       ifelse(ordered_samples %in% PRN4_BM, "BM", NA))

# Define annotation colors

group_colors <- c(LN = "grey", BM = "black")  # grey for LN, black for BM
celltype_colors <- c("Normal B-cell" = "#E7872B", "Burkitt Lymphoma cell" = "#3F78C1", "LN Tumour Bulk" = "#00008B")

# Define which samples are Normal B cell and which are bulk

normal_bcell_samples <- PRN4_rows$Sample_name[PRN4_rows$Myc_translocation_IGV == "No"]

PRN4_LN_Bulk <- c("PRN4GBDLBC72")

# Define CellType annotation

cell_type <- case_when(
  ordered_samples %in% normal_bcell_samples ~ "Normal B-cell",
  ordered_samples %in% PRN4_LN_Bulk        ~ "LN Tumour Bulk",
  TRUE                                      ~ "Burkitt Lymphoma cell"
)


# Build the annotation data.frame

annot_df <- data.frame(
  Location  = sample_group,   # BM / LN / NA
  Cell_Type = cell_type,      # factor we set earlier
  row.names = ordered_samples,
  check.names = FALSE
)

# Explicitly set factor levels so order() / arrange() know the hierarchy
annot_df$Cell_Type <- factor(
  annot_df$Cell_Type,
  levels = c("Normal B-cell",
             "Burkitt Lymphoma cell",
             "LN Tumour Bulk")   # overall order of the three blocks
)

# Within-BL block, BM first then LN

annot_df$Location <- factor(
  annot_df$Location,
  levels = c("BM", "LN")         
)

# Derive the new sample order
ordered_samples2 <- rownames(
  annot_df[ order(annot_df$Cell_Type, annot_df$Location), ]
)

#   Re-create the HeatmapAnnotation with the *reordered* annot_df
top_annot <- HeatmapAnnotation(
  df  = annot_df[ordered_samples2, ],   # keep rows aligned
  col = list(
    Location  = group_colors,
    Cell_Type = celltype_colors
  ),
  annotation_name_side = "left"
)

# Make initial oncoplot

pdf("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/PRN4_oncoprint_before_manual_check.pdf",
    width = 10, height = 6)

oncoPrint(
  effect_matrix_t[ , ordered_samples2, drop = FALSE],
  alter_fun            = alter_fun,
  col                  = col,
  column_title         = "PRN4 oncoplot",
  heatmap_legend_param = heatmap_legend_param,
  show_column_names    = TRUE,
  column_order         = ordered_samples2,
  top_annotation       = top_annot
)

dev.off()

# Make a table of driver mutations to check manually on IGV

drivers_to_check <- drivers_PRN4_shared %>%
  dplyr::select(-Sample) %>%
  distinct()

# Import csv file AFTER manual checking on IGV

effect_matrix_after_check <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/PRN4_manual_driver_check.csv")
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

missing_samples <- setdiff(PRN4_samples, colnames(effect_matrix_after_check_T))

# If there are missing samples, add them as new columns with "" values

new_cols <- matrix("", nrow = nrow(effect_matrix_after_check_T), ncol = length(missing_samples))
colnames(new_cols) <- missing_samples
rownames(new_cols) <- rownames(effect_matrix_after_check_T)
effect_matrix_after_check_T <- cbind(effect_matrix_after_check_T, new_cols)

# Make initial oncoplot

pdf("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/PRN4_oncoprint_after_manual_check.pdf",
    width = 10, height = 6)

oncoPrint(
  effect_matrix_after_check_T[ , ordered_samples2, drop = FALSE],
  alter_fun            = alter_fun,
  col                  = col,
  column_title         = "PRN4 oncoplot",
  heatmap_legend_param = heatmap_legend_param,
  show_column_names    = TRUE,
  column_order         = ordered_samples2,
  top_annotation       = top_annot
)

dev.off()
