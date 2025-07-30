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

# Load drivers found using the script: /hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus/P3G6

drivers_P3G6_batch_1 <- VariantAnnotation::readVcf("~//hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus/P3G6/P3G6_1.vep.effect.genes.somatic.vcf")
drivers_P3G6_batch_2 <- VariantAnnotation::readVcf("~//hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus/P3G6/P3G6_2.vep.effect.genes.somatic.vcf")

extract_genes_mutations <- function(vcf) {
  csq_field <- info(vcf)$CSQ
  
  # Get CSQ header format string
  csq_desc <- info(header(vcf))["CSQ", "Description"]
  csq_format <- sub(".*Format: ", "", csq_desc)
  csq_fields <- unlist(strsplit(csq_format, "\\|"))
  
  # Field indices
  gene_idx <- which(csq_fields %in% c("SYMBOL", "Gene", "Gene_Name"))[1]
  mut_idx  <- which(csq_fields %in% c("Consequence", "HGVSp", "HGVSc"))[1]  # Consequence first
  
  # Get variant IDs
  rr <- rowRanges(vcf)
  variant_ids <- paste0(
    seqnames(rr), ":", start(rr), "_",
    as.character(mcols(rr)$REF), ">",
    sapply(mcols(rr)$ALT, function(x) paste(as.character(x), collapse = ","))
  )
  
  # Flatten CSQ
  csq_flat <- unlist(csq_field)
  variant_index <- rep(seq_along(variant_ids), elementNROWS(csq_field))
  parsed_lines <- strsplit(csq_flat, "\\|")
  
  # Normalize length
  max_len <- max(lengths(parsed_lines))
  parsed_matrix <- t(sapply(parsed_lines, function(x) {
    length(x) <- max_len
    x
  }))
  
  # Output data frame
  df <- data.frame(
    variant = variant_ids[variant_index],
    gene = parsed_matrix[, gene_idx],
    mutation = parsed_matrix[, mut_idx],
    stringsAsFactors = FALSE
  )
  
  df <- df[df$gene != "", ]
  return(unique(df))
}
genes_mut_batch1 <- extract_genes_mutations(drivers_P3G6_batch_1)
genes_mut_batch1_dedup <- genes_mut_batch1[!duplicated(genes_mut_batch1$variant), ]
genes_mut_batch2 <- extract_genes_mutations(drivers_P3G6_batch_2)
genes_mut_batch2_dedup <- genes_mut_batch2[!duplicated(genes_mut_batch2$variant), ]
# Combine the two
genes_mut_merged <- rbind(genes_mut_batch1_dedup, genes_mut_batch2_dedup)
genes_mut_merged_unique <- genes_mut_merged[!duplicated(genes_mut_merged[, c("variant", "gene")]), ]
genes_mut_merged_unique_2 <- genes_mut_merged[!duplicated(genes_mut_merged$variant), ]

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
    Effect == "missense_variant"    ~ "Missense mutation",
    Effect == "missense_variant&splice_region_variant"    ~ "Missense mutation",
    Effect == "5_prime_UTR_variant"    ~ "UTR mutation",
    Effect == "3_prime_UTR_variant&NMD_transcript_variant"    ~ "UTR mutation",
    Effect == "stop_gained"         ~ "Nonsense mutation",
    Effect == "splice_donor_variant"   ~ "Splice site",
    Effect == "frameshift_variant"   ~ "Frameshift mutation",
    Effect == "downstream_gene_variant"    ~ "Up or Downstream gene mutation",
    Effect == "inframe_insertion"    ~ "In Frame Ins",
    Effect == "inframe_deletion&splice_region_variant"    ~ "In Frame Del",
    Effect == "upstream_gene_variant"    ~ "Up or Downstream gene mutation",
    Effect == "intron_variant"    ~ "Intron mutation",
    TRUE                            ~ Effect  # keep unchanged if not matched
  ))

# Make colour palette for later

col <- c(
  Missense_mutation              = "#008000",
  Intron_mutation                = "#7D3C98",
  In_Frame_Ins                   = "#8B0000",
  Frameshift_mutation            = "#1F77B4",
  Splice_site                    = "#FF8C00",
  Nonsense_mutation              = "#FF0000",
  UTR_mutation                   = "#FFD700",
  In_Frame_Del                   = "#8B4513",
  Up_or_Downstream_gene_mutation    = "#20B2AA"
)

# Make a table of driver mutations to check manually on IGV

drivers_to_check <- drivers_P3G6_shared %>% dplyr::select(-Sample) %>% distinct()

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

# Provide sample info

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
P3G6_rows <- input_df %>% filter(Novogene_ID == "P3G6") %>% arrange(factor(Myc_translocation_IGV, levels = c("no", "yes")))   
P3G6_rows <- P3G6_rows %>% filter(!grepl("MSCBULK", Sample_name))
P3G6_samples <- P3G6_rows$Sample_name

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

normal_bcell_samples <- P3G6_rows$Sample_name[P3G6_rows$Myc_translocation_IGV == "No"]

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

# Create boxes

draw_box <- function(x, y, w, h, fill)
  grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
            gp = gpar(fill = fill, col = NA))

variant_classes <- names(col)

alter_fun <- c(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  setNames(lapply(variant_classes, \(v)
                  \(x, y, w, h) draw_box(x, y, w, h, col[v])), variant_classes)
)

# Create the top annotation with both Location and CellType

top_annot <- HeatmapAnnotation(
  Location = sample_group,
  Sample = cell_type,
  col = list(
    Location = group_colors,
    Sample = celltype_colors
  ),
  annotation_name_side = "left"
)

column_title = "P3G6 oncoplot"

# Variant labels

heatmap_legend_param <- list(
  title  = "Alterations",
  at     = names(col),    
  labels = gsub("_", " ", names(col))
)

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

