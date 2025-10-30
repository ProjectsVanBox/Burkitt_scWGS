################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to plot oncoplot for PIA9
# Author: Alexander Steemers
# Date: August 2025
################################################################################

# Load libraries

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(readxl)
library(purrr)
library(VariantAnnotation)
library(circlize) 
library(grid)

# Set working directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/")

# Find blacklisted samples

low_callable_df<- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv')
below_curve_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df     <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df    <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')

blacklist <- unique(c(below_curve_df$Sample_name,
                      low_callable_df$Sample_name,
                      bad_baf_df$Sample_name,
                      fail_vaf_df$samplename))

# Get known BL driver genes

BL_drivers <- read.delim(
  "~/hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus/BL_Intogen_Panea_Blood2019_Nicole_Blood_2023.txt",
  header = FALSE, stringsAsFactors = FALSE
)[,1]

# Get known cancer drivers

Cancer_drivers <- read.delim(
  "~/hpc/pmc_vanboxtel/resources/hotspot_genes/cancerGeneList_oncokb_06262025_genenames.txt",
  header = FALSE, stringsAsFactors = FALSE
)[,1] 

# Combine

all_known_drivers <- unique(c(BL_drivers, Cancer_drivers))

# Function to extract driver mutations

extract_genes_mutations <- function(vcf) {
  # --- 1) Parse CSQ field layout from header ---
  csq_field  <- VariantAnnotation::info(vcf)$CSQ
  csq_desc   <- as.character(VariantAnnotation::info(header(vcf))["CSQ", "Description"])
  csq_format <- sub(".*Format:\\s*", "", csq_desc)
  csq_fields <- strsplit(csq_format, "\\|")[[1]]
  
  # indices of useful fields
  gene_idx <- which(csq_fields %in% c("SYMBOL","Gene","Gene_Name"))[1]
  conseq_idx <- which(csq_fields == "Consequence")[1]
  hgvsp_idx  <- which(csq_fields == "HGVSp")[1]
  hgvsc_idx  <- which(csq_fields == "HGVSc")[1]
  
  # --- 2) Variant identifiers ---
  rr   <- rowRanges(vcf)
  REF  <- as.character(mcols(rr)$REF)
  ALTl <- mcols(rr)$ALT
  ALT  <- vapply(ALTl, function(x) paste(as.character(x), collapse=","), "")
  
  chrom <- as.character(seqnames(rr))
  pos   <- start(rr)
  variant_ids <- paste0(chrom, ":", pos, "_", REF, ">", ALT)
  
  # --- 3) Flatten CSQ annotations to matrix ---
  csq_flat   <- unlist(csq_field, use.names = FALSE)
  var_index  <- rep(seq_along(variant_ids), lengths(csq_field))
  parsed     <- strsplit(csq_flat, "\\|")
  max_len    <- max(vapply(parsed, length, 0L))
  parsed_mat <- t(vapply(parsed, function(x){ length(x) <- max_len; x }, character(max_len)))
  
  ann_df <- data.frame(
    variant   = variant_ids[var_index],
    chrom     = chrom[var_index],
    pos       = pos[var_index],
    ref       = REF[var_index],
    alt       = ALT[var_index],
    gene      = if (!is.na(gene_idx)) parsed_mat[, gene_idx] else NA,
    consequence = if (!is.na(conseq_idx)) parsed_mat[, conseq_idx] else NA,
    hgvsp     = if (!is.na(hgvsp_idx)) parsed_mat[, hgvsp_idx] else NA,
    hgvsc     = if (!is.na(hgvsc_idx)) parsed_mat[, hgvsc_idx] else NA,
    stringsAsFactors = FALSE
  )
  ann_df <- ann_df[nzchar(ann_df$gene), , drop = FALSE]
  
  # --- 4) Pull AD and compute ref/alt counts + VAF ---
  AD <- geno(vcf)$AD
  if (is.null(AD)) stop("FORMAT/AD not present in VCF; cannot compute VAF.")
  
  samples <- colnames(AD)
  counts_long <- do.call(rbind, lapply(seq_len(nrow(AD)), function(i) {
    cell_list <- as.list(AD[i, , drop=TRUE])
    do.call(rbind, lapply(seq_along(cell_list), function(j) {
      cnt <- as.integer(cell_list[[j]])
      ref_cnt <- if (length(cnt) >= 1) cnt[1] else NA
      alt_cnt <- if (length(cnt) >= 2) sum(cnt[-1]) else 0L
      dp <- ref_cnt + alt_cnt
      data.frame(
        variant   = variant_ids[i],
        sample    = samples[j],
        ref_count = ref_cnt,
        alt_count = alt_cnt,
        dp        = dp,
        vaf       = if (!is.na(dp) && dp > 0) alt_cnt/dp else NA,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # --- 5) Join ---
  out <- merge(ann_df, counts_long, by="variant", all.x=TRUE)
  out <- unique(out[ , c("variant","chrom","pos","ref","alt",
                         "gene","consequence","hgvsp","hgvsc",
                         "sample","ref_count","alt_count","dp","vaf") ])
  out
}

# folder containing the PRN4 batch VCFs

vcf_dir <- "~/hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus/PIA9/Single_cell_drivers"

# find all batches (handles both underscore and dot separators)

vcf_files <- list.files(
  vcf_dir,
  pattern = ".vep\\.effect_3\\.somatic\\.pass\\.vcf$",
  full.names = TRUE
)

if (length(vcf_files) == 0) stop("No PIA9 batch VCFs found in: ", vcf_dir)

# read all VCFs, extract tables, and combine

genes_mut_merged <- vcf_files %>%
  map(~ VariantAnnotation::readVcf(.x)) %>%
  map(extract_genes_mutations) %>%
  bind_rows()

# Remove duplicate rows

genes_mut_merged_nodup <- unique(genes_mut_merged)

# If a particular variant has multiple mutation annotations keep only the first annotation

#genes_mut_merged_filterd_1 <- genes_mut_merged_nodup %>% distinct(variant, sample, .keep_all = TRUE)

# Remove rows with samples which have been blacklisted

genes_mut_merged_filterd_2 <- genes_mut_merged_nodup %>%
  filter(!sample %in% blacklist)

# Filter out rows which do not have genes in the BL list of known genes

genes_mut_merged_filterd_3 <- genes_mut_merged_filterd_2 %>% filter(gene %in% BL_drivers)

# Simplify mutation annotation

genes_mut_merged_filterd_4 <- genes_mut_merged_filterd_3 %>%
  mutate(consequence = stringr::str_squish(stringr::str_replace(consequence, "&.*$", "")))


# Remove rows with VAF 0 or N/A

genes_mut_merged_filterd_5 <- genes_mut_merged_filterd_4 %>%
  filter(!is.na(vaf), vaf != 0)

# Remove rows with protein-altering variants (i.e. discard synonymous, intron variants, up/downstream variants etc.)

driver_terms <- c(
  "missense_variant",
  "stop_gained",
  "stop_lost",
  "start_lost",
  "frameshift_variant",
  "inframe_insertion",
  "inframe_deletion",
  "protein_altering_variant",
  "coding_sequence_variant"
)

genes_mut_merged_filterd_6 <- subset(genes_mut_merged_filterd_5, grepl(paste(driver_terms, collapse="|"), consequence))

# Provide sample info

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
PIA9_rows <- input_df %>% filter(Novogene_ID == "PIA9") %>% arrange(factor(Myc_translocation_IGV, levels = c("no", "yes")))   
normal_bcell_samples <- PIA9_rows$Sample_name[PIA9_rows$Myc_translocation_IGV == "No"]

# Remove drivers if they are found in normal b cells

genes_mut_merged_filterd_7 <- genes_mut_merged_filterd_6 %>%
  group_by(variant) %>%
  filter(sum(sample %in% normal_bcell_samples) < 1) %>%
  ungroup()

genes_mut_merged_filterd_8 <- genes_mut_merged_filterd_7 %>%
  mutate(driver = paste0(gene, "_", variant))

genes_mut_merged_filterd_9 <- genes_mut_merged_filterd_8 %>%
  dplyr::group_by(driver, sample) %>%
  dplyr::slice(1) %>%   # keep first row per group
  dplyr::ungroup()

genes_mut_merged_filterd_10 <- unique(genes_mut_merged_filterd_9)

# Remove variants that only occur once

genes_mut_merged_filterd_11 <- genes_mut_merged_filterd_10 %>%
  add_count(variant) %>%
  filter(n > 1) %>%
  dplyr::select(-n) 

# Prepare oncoplot matrix 

oncoplot_matrix <- genes_mut_merged_filterd_11 %>%
  dplyr::select(driver, sample, consequence) %>%      # keep only relevant columns
  pivot_wider(
    names_from = sample,
    values_from = consequence,
    values_fill = list(mutation = NA)        # fill missing with NA
  )


# Desired column order (from CellPhy tree)

# Desired column order (from CellPhy tree)

PIA9_samples <- c(
  "PIA9GTDBBC75","PIA9GTDBBC64","PIA9GTDBBC52","PIA9GTDBBC54","PIA9GTDBBC57",
  "PIA9GTDBBC73","PIA9GTDBBC58","PIA9GTDBBC60","PIA9GTDBBC67","PIA9GTDBBC59",
  "PIA9GTDBBC63","PIA9GTDBBC72","PIA9GTDBBC61","PIA9GTDBBC77","PIA9GTDBBC68",
  "PIA9GTDBBC55","PIA9GTDBBC65","PIA9GTDBBC53","PIA9GTDBBC66","PIA9GTDBBC56",
  "PIA9GTDBBC74","PIA9GTDBBC71","PIA9GTDBBC76",
  "PIA9GTDABC37","PIA9GTDABC7","PIA9GTDABC49","PIA9GTDABC21","PIA9GTDABC33",
  "PIA9GTDABC19","PIA9GTDABC9","PIA9GTDABC26","PIA9GTDABC44","PIA9GTDABC5",
  "PIA9GTDABC6","PIA9GTDABC25","PIA9GTDABC48","PIA9GTDABC18","PIA9GTDABC40",
  "PIA9GTDABC34","PIA9GTDABC50","PIA9GTDABC42","PIA9GTDABC46","PIA9GTDABC47",
  "PIA9GTDABC43","PIA9GTDABC20","PIA9GTDABC36","PIA9GTDABC24","PIA9GTDABC39",
  "PIA9GTDABC15","PIA9GTDABC51","PIA9GTDABC23","PIA9GTDABC17","PIA9GTDABC45",
  "PIA9GTDABC41","PIA9GTDABC27","PIA9GTDABC35","PIA9GTDABC22"
)


# Ensure all PIA9_samples exist as columns (add NA if missing)

for (s in PIA9_samples) {
  if (!s %in% colnames(oncoplot_matrix)) {
    oncoplot_matrix[[s]] <- NA
  }
}

# Reorder columns: driver first, then PIA9_samples
oncoplot_matrix <- oncoplot_matrix %>%
  dplyr::select(driver, all_of(PIA9_samples))

# export driver list 

drivers_by_sample_PIA9 <- oncoplot_matrix %>%
  pivot_longer(
    cols = -driver,
    names_to = "sample",
    values_to = "val",
    values_transform = list(val = as.character)  # coerce mixed types safely
  ) %>%
  # keep only calls that indicate presence (not NA/empty/"0"/"false")
  filter(!is.na(val),
         str_trim(val) != "",
         !val %in% c("NA", "0"),
         tolower(val) != "false") %>%
  dplyr::select(sample, driver) %>%
  distinct() %>%                 # avoid duplicates if any
  arrange(sample, driver)

write_xlsx(drivers_by_sample_PIA9, "Data/drivers_by_sample_PIA9.xlsx")

# Build the effect matrix from oncoplot_matrix

effect_df <- oncoplot_matrix

# Keep driver column during wrangling
effect_df <- effect_df %>% relocate(driver)

# Right before converting to matrix:

effect_matrix_after_check_T <- as.matrix(effect_df[,-1])
rownames(effect_matrix_after_check_T) <- effect_df$driver

# Replace NAs by empty strings (ComplexHeatmap uses "" as absence)
effect_matrix_after_check_T[is.na(effect_matrix_after_check_T)] <- ""

# Column order from sample sheet (and keep only columns that exist)

ordered_samples <- colnames(effect_matrix_after_check_T)

# Make a filtered sample list that (a) isn’t blacklisted and (b) exists in the matrix
PIA9_samples_filtered <- setdiff(PIA9_samples, blacklist)
PIA9_samples_filtered <- intersect(PIA9_samples_filtered, colnames(effect_matrix_after_check_T))

# Cell type colors and vector aligned to ordered_samples
celltype_colors <- c("Normal B-cell" = "#E7872B",
                     "Burkitt Lymphoma cell" = "#3F78C1",
                     "Bulk Tumour" = "#00008B")

# Categorize samples in PIA9_rows
PIA9_rows <- PIA9_rows %>%
  mutate(
    Category = case_when(
      Sample_name %in% normal_bcell_samples ~ "Normal B-cell",
      grepl("BULK", Sample_name, ignore.case = TRUE) ~ "Bulk Tumour",
      TRUE ~ "Burkitt Lymphoma cell"
    ),
    Category = factor(Category,
                      levels = c("Normal B-cell", "Burkitt Lymphoma cell", "Bulk Tumour"))
  )

# Make a named lookup and align it to ordered_samples
sample_to_cat <- setNames(as.character(PIA9_rows$Category), PIA9_rows$Sample_name)
cell_type <- unname(sample_to_cat[ordered_samples])

# Variant class colors (ensure these strings match your 'mutation' values)
col <- c(
  missense_variant    = "#008000",   # green
  intron_variant      = "#7D3C98",   # purple
  inframe_insertion   = "#8B0000",   # dark red
  frameshift_variant  = "#1F77B4",   # blue
  protein_altering_variant         = "#FF8C00",   # orange
  stop_gained         = "#FF0000",   # red
  "3_prime_UTR_variant&NMD_transcript_variant" = "#FFD700", # gold
  inframe_deletion    = "#8B4513",   # brown
  Up_or_Downstream_gene_mutation = "#20B2AA",   # teal
  start_lost = "#A9A9A9",   # dark gray
  downstream_gene_variant            = "#40E0D0",   # turquoise
  "frameshift_variant&splice_region_variant" = "#6495ED", # cornflower blue
  "intron_variant&NMD_transcript_variant" = "#BA55D3", # medium orchid
  upstream_gene_variant              = "#00CED1",   # dark turquoise
  "5_prime_UTR_variant"              = "#DAA520",   # goldenrod
  "inframe_insertion&stop_retained_variant" = "#CD5C5C", # indian red
  "splice_region_variant&splice_polypyrimidine_tract_variant&intron_variant" = "#FF69B4", # hot pink
  "frameshift_variant&stop_lost"     = "#4682B4"   , # steel blue
  "3_prime_UTR_variant" = "black" # black
)


# Alteration drawing functions
draw_box <- function(x, y, w, h, fill)
  grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"), gp = gpar(fill = fill, col = NA))

variant_classes <- names(col)
alter_fun <- c(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  setNames(lapply(variant_classes, \(v) \(x, y, w, h) draw_box(x, y, w, h, col[v])), variant_classes)
)

# Top annotation aligned to ordered_samples
top_annot <- HeatmapAnnotation(
  Sample = cell_type,
  col = list(
    Sample = celltype_colors
  ),
  annotation_name_side = "left"
)

# Legend labels
heatmap_legend_param <- list(
  title  = "Alterations",
  at     = names(col),
  labels = gsub("_", " ", names(col))
)

column_title <- "PIA9 oncoplot"

# Plot
pdf("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/Figures/PIA9_oncoprint_burkitt_genes.pdf",
    width = 10, height = 6)

oncoPrint(
  effect_matrix_after_check_T[, ordered_samples, drop = FALSE],
  alter_fun = alter_fun,
  col = col,
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param,
  show_column_names = TRUE,
  column_order = ordered_samples,
  top_annotation = top_annot,
  right_annotation = NULL,
  remove_empty_rows = TRUE,
  remove_empty_columns = FALSE,
  column_names_gp = gpar(fontsize = 10),  
  row_names_gp   = gpar(fontsize = 10)
)

dev.off()

# Filter out rows which do not have genes in the BL list of known genes

genes_mut_merged_filterd_3 <- genes_mut_merged_filterd_2 %>% filter(gene %in% all_known_drivers)

# Simplify mutation annotation

genes_mut_merged_filterd_4 <- genes_mut_merged_filterd_3 %>%
  mutate(consequence = stringr::str_squish(stringr::str_replace(consequence, "&.*$", "")))


# Remove rows with VAF 0 or N/A

genes_mut_merged_filterd_5 <- genes_mut_merged_filterd_4 %>%
  filter(!is.na(vaf), vaf != 0)

# Remove rows with protein-altering variants (i.e. discard synonymous, intron variants, up/downstream variants etc.)

driver_terms <- c(
  "missense_variant",
  "stop_gained",
  "stop_lost",
  "start_lost",
  "frameshift_variant",
  "inframe_insertion",
  "inframe_deletion",
  "protein_altering_variant",
  "coding_sequence_variant"
)

genes_mut_merged_filterd_6 <- subset(genes_mut_merged_filterd_5, grepl(paste(driver_terms, collapse="|"), consequence))

# Provide sample info

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
PIA9_rows <- input_df %>% filter(Novogene_ID == "PIA9") %>% arrange(factor(Myc_translocation_IGV, levels = c("no", "yes")))   
normal_bcell_samples <- PIA9_rows$Sample_name[PIA9_rows$Myc_translocation_IGV == "No"]

# Remove drivers if they are found in normal b cells

genes_mut_merged_filterd_7 <- genes_mut_merged_filterd_6 %>%
  group_by(variant) %>%
  filter(sum(sample %in% normal_bcell_samples) < 1) %>%
  ungroup()

genes_mut_merged_filterd_8 <- genes_mut_merged_filterd_7 %>%
  mutate(driver = paste0(gene, "_", variant))

genes_mut_merged_filterd_9 <- genes_mut_merged_filterd_8 %>%
  dplyr::group_by(driver, sample) %>%
  dplyr::slice(1) %>%   # keep first row per group
  dplyr::ungroup()

genes_mut_merged_filterd_10 <- unique(genes_mut_merged_filterd_9)

# Remove variants that only occur once

genes_mut_merged_filterd_11 <- genes_mut_merged_filterd_10 %>%
  add_count(variant) %>%
  filter(n > 1) %>%
  dplyr::select(-n) 

# Prepare oncoplot matrix 

oncoplot_matrix <- genes_mut_merged_filterd_11 %>%
  dplyr::select(driver, sample, consequence) %>%      # keep only relevant columns
  pivot_wider(
    names_from = sample,
    values_from = consequence,
    values_fill = list(mutation = NA)        # fill missing with NA
  )

# Ensure all PIA9_samples exist as columns (add NA if missing)

for (s in PIA9_samples) {
  if (!s %in% colnames(oncoplot_matrix)) {
    oncoplot_matrix[[s]] <- NA
  }
}

# Reorder columns: driver first, then PIA9_samples
oncoplot_matrix <- oncoplot_matrix %>%
  dplyr::select(driver, all_of(PIA9_samples))

# Build the effect matrix from oncoplot_matrix

effect_df <- oncoplot_matrix

# Keep driver column during wrangling
effect_df <- effect_df %>% relocate(driver)

# Right before converting to matrix:

effect_matrix_after_check_T <- as.matrix(effect_df[,-1])
rownames(effect_matrix_after_check_T) <- effect_df$driver

# Replace NAs by empty strings (ComplexHeatmap uses "" as absence)
effect_matrix_after_check_T[is.na(effect_matrix_after_check_T)] <- ""

# Column order from sample sheet (and keep only columns that exist)

ordered_samples <- colnames(effect_matrix_after_check_T)

# Make a filtered sample list that (a) isn’t blacklisted and (b) exists in the matrix
PIA9_samples_filtered <- setdiff(PIA9_samples, blacklist)
PIA9_samples_filtered <- intersect(PIA9_samples_filtered, colnames(effect_matrix_after_check_T))

# Cell type colors and vector aligned to ordered_samples
celltype_colors <- c("Normal B-cell" = "#E7872B",
                     "Burkitt Lymphoma cell" = "#3F78C1",
                     "Bulk Tumour" = "#00008B")

# Categorize samples in PIA9_rows
PIA9_rows <- PIA9_rows %>%
  mutate(
    Category = case_when(
      Sample_name %in% normal_bcell_samples ~ "Normal B-cell",
      grepl("BULK", Sample_name, ignore.case = TRUE) ~ "Bulk Tumour",
      TRUE ~ "Burkitt Lymphoma cell"
    ),
    Category = factor(Category,
                      levels = c("Normal B-cell", "Burkitt Lymphoma cell", "Bulk Tumour"))
  )

# Make a named lookup and align it to ordered_samples
sample_to_cat <- setNames(as.character(PIA9_rows$Category), PIA9_rows$Sample_name)
cell_type <- unname(sample_to_cat[ordered_samples])

# Variant class colors (ensure these strings match your 'mutation' values)
col <- c(
  missense_variant    = "#008000",   # green
  intron_variant      = "#7D3C98",   # purple
  inframe_insertion   = "#8B0000",   # dark red
  frameshift_variant  = "#1F77B4",   # blue
  protein_altering_variant         = "#FF8C00",   # orange
  stop_gained         = "#FF0000",   # red
  "3_prime_UTR_variant&NMD_transcript_variant" = "#FFD700", # gold
  inframe_deletion    = "#8B4513",   # brown
  Up_or_Downstream_gene_mutation = "#20B2AA",   # teal
  start_lost = "#A9A9A9",   # dark gray
  downstream_gene_variant            = "#40E0D0",   # turquoise
  "frameshift_variant&splice_region_variant" = "#6495ED", # cornflower blue
  "intron_variant&NMD_transcript_variant" = "#BA55D3", # medium orchid
  upstream_gene_variant              = "#00CED1",   # dark turquoise
  "5_prime_UTR_variant"              = "#DAA520",   # goldenrod
  "inframe_insertion&stop_retained_variant" = "#CD5C5C", # indian red
  "splice_region_variant&splice_polypyrimidine_tract_variant&intron_variant" = "#FF69B4", # hot pink
  "frameshift_variant&stop_lost"     = "#4682B4"   , # steel blue
  "3_prime_UTR_variant" = "black" # black
)


# Alteration drawing functions
draw_box <- function(x, y, w, h, fill)
  grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"), gp = gpar(fill = fill, col = NA))

variant_classes <- names(col)
alter_fun <- c(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(2, "pt"), h - unit(2, "pt"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  setNames(lapply(variant_classes, \(v) \(x, y, w, h) draw_box(x, y, w, h, col[v])), variant_classes)
)

# Top annotation aligned to ordered_samples
top_annot <- HeatmapAnnotation(
  Sample = cell_type,
  col = list(
    Sample = celltype_colors
  ),
  annotation_name_side = "left"
)

# Legend labels
heatmap_legend_param <- list(
  title  = "Alterations",
  at     = names(col),
  labels = gsub("_", " ", names(col))
)

column_title <- "PIA9 oncoplot"

# Plot
pdf("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/Figures/PIA9_oncoprint_pan_cancer_genes.pdf",
    width = 10, height = 6)

oncoPrint(
  effect_matrix_after_check_T[, ordered_samples, drop = FALSE],
  alter_fun = alter_fun,
  col = col,
  column_title = column_title,
  heatmap_legend_param = heatmap_legend_param,
  show_column_names = TRUE,
  column_order = ordered_samples,
  top_annotation = top_annot,
  right_annotation = NULL,
  remove_empty_rows = TRUE,
  remove_empty_columns = FALSE,
  column_names_gp = gpar(fontsize = 6),  
  row_names_gp   = gpar(fontsize = 6)
)

dev.off()