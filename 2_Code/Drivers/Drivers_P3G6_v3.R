################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to plot oncoplot for P3G6
# Author: Alexander Steemers
# Date: August 2025
################################################################################

# Load libraries

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(readxl)
library(purrr)
library(ComplexHeatmap)
library(circlize) 
library(grid)
library(VariantAnnotation)

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
  "~/hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus/burkitt_drivers.txt",
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
  csq_field  <- info(vcf)$CSQ
  csq_desc   <- as.character(info(header(vcf))["CSQ", "Description"])
  csq_format <- sub(".*Format:\\s*", "", csq_desc)
  csq_fields <- strsplit(csq_format, "\\|")[[1]]
  
  # indices of fields we care about
  gene_idx <- which(csq_fields %in% c("SYMBOL", "Gene", "Gene_Name"))[1]
  mut_idx  <- which(csq_fields %in% c("Consequence", "HGVSp", "HGVSc"))[1]  # prefer Consequence
  
  # --- 2) Variant identifiers (CHROM:POS_REF>ALT) + REF/ALT strings ---
  rr   <- rowRanges(vcf)
  REF  <- as.character(mcols(rr)$REF)
  ALTl <- mcols(rr)$ALT                       # DNAStringSetList
  ALT  <- vapply(ALTl, function(x) paste(as.character(x), collapse=","), "")
  
  chrom <- as.character(seqnames(rr))
  pos   <- start(rr)
  variant_ids <- paste0(chrom, ":", pos, "_", REF, ">", ALT)
  
  # --- 3) Flatten CSQ annotations to a matrix ---
  csq_flat   <- unlist(csq_field, use.names = FALSE)
  var_index  <- rep(seq_along(variant_ids), lengths(csq_field))  # which variant each CSQ row belongs to
  parsed     <- strsplit(csq_flat, "\\|")
  max_len    <- max(vapply(parsed, length, 0L))
  parsed_mat <- t(vapply(parsed, function(x){ length(x) <- max_len; x }, character(max_len)))
  
  gene_col <- if (!is.na(gene_idx)) parsed_mat[, gene_idx] else rep(NA_character_, nrow(parsed_mat))
  mut_col  <- if (!is.na(mut_idx))  parsed_mat[, mut_idx]  else rep(NA_character_, nrow(parsed_mat))
  
  ann_df <- data.frame(
    variant = variant_ids[var_index],
    chrom   = chrom[var_index],
    pos     = pos[var_index],
    ref     = REF[var_index],
    alt     = ALT[var_index],
    gene    = gene_col,
    mutation= mut_col,
    stringsAsFactors = FALSE
  )
  ann_df <- ann_df[nzchar(ann_df$gene), , drop = FALSE]
  
  # --- 4) Pull AD and compute per-sample ref/alt counts + VAF ---
  # AD is a matrix (variants x samples), each cell = integer vector c(REF, ALT1, ALT2, ...)
  AD <- geno(vcf)$AD
  if (is.null(AD)) stop("FORMAT/AD not present in VCF; cannot compute VAF.")
  
  samples <- colnames(AD)
  # build long table of counts for all variants x samples
  counts_long <- do.call(rbind, lapply(seq_len(nrow(AD)), function(i) {
    cell_list <- as.list(AD[i, , drop = TRUE])
    do.call(rbind, lapply(seq_along(cell_list), function(j) {
      cnt <- as.integer(cell_list[[j]])
      ref_cnt <- if (length(cnt) >= 1) cnt[1] else NA_integer_
      alt_cnt <- if (length(cnt) >= 2) sum(cnt[-1]) else 0L     # sum all ALT alleles if multi-allelic
      dp      <- ref_cnt + alt_cnt
      data.frame(
        variant = variant_ids[i],
        sample  = samples[j],
        ref_count = ref_cnt,
        alt_count = alt_cnt,
        dp        = dp,
        vaf       = if (!is.na(dp) && dp > 0) alt_cnt / dp else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  # --- 5) Join CSQ annotations with counts (many-to-many if multiple CSQs per variant) ---
  out <- merge(ann_df, counts_long, by = "variant", all.x = TRUE)
  
  # unique rows, keep useful ordering
  out <- unique(out[ , c("variant","chrom","pos","ref","alt","gene","mutation",
                         "sample","ref_count","alt_count","dp","vaf") ])
  
  out
}

# folder containing the P3G6 batch VCFs

vcf_dir <- "~/hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus/P3G6/Single_cell_drivers"

# find all batches (handles both underscore and dot separators)

vcf_files <- list.files(
  vcf_dir,
  #pattern = "^P3G6[._][0-9]+\\.vep\\.effect\\.genes\\.somatic\\.pass\\.vcf$",
  #pattern = "^P3G6[._][0-9]+\\.vep\\.effect_2\\.genes_2\\.somatic_2\\.pass_2\\.vcf$",
  pattern = ".vep\\.effect_3\\.somatic\\.pass\\.vcf$",
  full.names = TRUE
)

if (length(vcf_files) == 0) stop("No P3G6 batch VCFs found in: ", vcf_dir)

# read all VCFs, extract tables, and combine

genes_mut_merged <- vcf_files %>%
  map(~ VariantAnnotation::readVcf(.x)) %>%
  map(extract_genes_mutations) %>%
  bind_rows()

# Remove duplicate rows

genes_mut_merged_nodup <- unique(genes_mut_merged)

# If a particular variant has multiple mutation annotations keep only the first annotation

genes_mut_merged_filterd_1 <- genes_mut_merged_nodup %>% distinct(variant, sample, .keep_all = TRUE)

# Remove rows with samples which have been blacklisted

genes_mut_merged_filterd_2 <- genes_mut_merged_filterd_1 %>%
  filter(!sample %in% blacklist)

# Filter out rows which do not have genes in the BL list of known genes

genes_mut_merged_filterd_3 <- genes_mut_merged_filterd_2 #%>%
  #filter(gene %in% BL_drivers)

# Simplify mutation annotation

genes_mut_merged_filterd_4 <- genes_mut_merged_filterd_3 %>%
  mutate(mutation = stringr::str_squish(stringr::str_replace(mutation, "&.*$", "")))


# Remove rows with VAF 0 or N/A

genes_mut_merged_filterd_5 <- genes_mut_merged_filterd_4 %>%
  filter(!is.na(vaf), vaf != 0)

# Remove rows with protein-altering variants (i.e. discard synonymous, intron variants, up/downstream variants etc.)

driver_terms <- c("missense_variant",
                  "frameshift_variant",
                  "inframe_deletion",
                  "stop_gained",
                  "inframe_insertion",
                  "protein_altering_variant")

genes_mut_merged_filterd_6 <- subset(genes_mut_merged_filterd_5, grepl(paste(driver_terms, collapse="|"), mutation))

# Remove variants that only occur once

genes_mut_merged_filterd_7 <- genes_mut_merged_filterd_6 %>%
  add_count(variant) %>%
  filter(n > 1) %>%
  dplyr::select(-n) 

# Provide sample info

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
P3G6_rows <- input_df %>% filter(Novogene_ID == "P3G6") %>% arrange(factor(Myc_translocation_IGV, levels = c("no", "yes")))   
normal_bcell_samples <- P3G6_rows$Sample_name[P3G6_rows$Myc_translocation_IGV == "No"]

# Remove drivers if they are found in normal b cells

filtered_df_final <- genes_mut_merged_filterd_7 %>%
  group_by(variant) %>%
  filter(sum(sample %in% normal_bcell_samples) < 1) %>%
  ungroup()

filtered_df_final <- filtered_df_final %>%
  mutate(driver = paste0(gene, "_", variant))

# Prepare oncoplot matrix 

oncoplot_matrix <- filtered_df_final %>%
  dplyr::select(driver, sample, mutation) %>%      # keep only relevant columns
  pivot_wider(
    names_from = sample,
    values_from = mutation,
    values_fill = list(mutation = NA)        # fill missing with NA
  )


# Desired column order (from CellPhy tree)
P3G6_samples <- c(
  "P3G6GPDABC31", "PB11197-BLASC-BCELLP2F4", "PB11197-BLASC-BCELLP2D4", 
  "PB11197-BLASC-BCELLP2C4", "PB11197-BLASC-BCELLP2B4", "PB11197-BLASC-BCELLP2E4", 
  "P3G6GPDABC28", "PB11197-BLASC-BCELLP1P3", "PB11197-BLASC-BCELLP1C4", 
  "PB11197-BLASC-BCELLP1L3", "PB11197-BLASC-BCELLP1J3", "PB11197-BLASC-BCELLP1K4", 
  "PB11197-BLASC-BCELLP1O3", "PB11197-BLASC-BCELLP1I4", "PB11197-BLASC-BCELLP1B4", 
  "P3G6GPDABC26"
)

# Ensure all P3G6_samples exist as columns (add NA if missing)

for (s in P3G6_samples) {
  if (!s %in% colnames(oncoplot_matrix)) {
    oncoplot_matrix[[s]] <- NA
  }
}

# Reorder columns: driver first, then P3G6_samples
oncoplot_matrix <- oncoplot_matrix %>%
  dplyr::select(driver, all_of(P3G6_samples))

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
P3G6_samples_filtered <- setdiff(P3G6_samples, blacklist)
P3G6_samples_filtered <- intersect(P3G6_samples_filtered, colnames(effect_matrix_after_check_T))

# Group colors and group vector (aligned to ordered_samples)
group_colors <- c(Ascites = "#000000")
sample_group <- ifelse(ordered_samples %in% P3G6_samples_filtered, "Ascites", NA)

# Cell type colors and vector aligned to ordered_samples
celltype_colors <- c("Normal B-cell" = "#E7872B",
                     "Burkitt Lymphoma cell" = "#3F78C1",
                     "Bulk Tumour" = "#00008B")

# Categorize samples in P3G6_rows
P3G6_rows <- P3G6_rows %>%
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
sample_to_cat <- setNames(as.character(P3G6_rows$Category), P3G6_rows$Sample_name)
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
  non_coding_transcript_exon_variant = "#A9A9A9",   # dark gray
  downstream_gene_variant            = "#40E0D0",   # turquoise
  "frameshift_variant&splice_region_variant" = "#6495ED", # cornflower blue
  "intron_variant&NMD_transcript_variant" = "#BA55D3", # medium orchid
  upstream_gene_variant              = "#00CED1",   # dark turquoise
  "5_prime_UTR_variant"              = "#DAA520",   # goldenrod
  "inframe_insertion&stop_retained_variant" = "#CD5C5C", # indian red
  "splice_region_variant&splice_polypyrimidine_tract_variant&intron_variant" = "#FF69B4", # hot pink
  "frameshift_variant&stop_lost"     = "#4682B4"    # steel blue
)

# Annotate from which list each driver comes from 

driver_to_gene <- setNames(filtered_df_final$gene, filtered_df_final$driver)
row_genes <- unname(driver_to_gene[rownames(effect_matrix_after_check_T)])

is_bl     <- row_genes %in% BL_drivers
is_cancer <- row_genes %in% Cancer_drivers

driver_class <- ifelse(is_bl & is_cancer, "Both",
                       ifelse(is_bl, "Burkitt_list_only",
                              ifelse(is_cancer, "Pan_cancer_list_only", "Neither")))

driver_class_cols <- c(
  Burkitt_list_only = "#FF0000",
  Pan_cancer_list_only  = "#0000FF",
  Both         = "#00A000",
  Neither      = "#CCCCCC"
)

# Row annotation (left side)
row_annot <- rowAnnotation(
  DriverClass = driver_class,
  col = list(DriverClass = driver_class_cols),
  annotation_name_side = "top"
)

row_label_cols <- driver_class_cols[driver_class]

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
  Location = sample_group,
  Sample = cell_type,
  col = list(
    Location = group_colors,
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

column_title <- "P3G6 oncoplot"

# Plot
pdf("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/Figures/P3G6_oncoprint_20250926.pdf",
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
  left_annotation = row_annot, 
  right_annotation = NULL,
  remove_empty_rows = TRUE,
  remove_empty_columns = FALSE
)

dev.off()


samples_in_clade_wo_driver <- c(
  "P3G6GPDABC26",
  "PB11197-BLASC-BCELLP1B",
  "PB11197-BLASC-BCELLP1I4",
  "PB11197-BLASC-BCELLP1O3",
  "PB11197-BLASC-BCELLP1K4",
  "PB11197-BLASC-BCELLP1J3"
)
# Filter drivers to the clade-only (>=5)

# Presence/absence matrix (TRUE if any alteration recorded)

presence <- effect_matrix_after_check_T != ""

# Define target-vs-other sample sets

clade_samples <- samples_in_clade_wo_driver
in_clade      <- colnames(effect_matrix_after_check_T) %in% clade_samples
in_other      <- !in_clade

# Counts per driver (row)

count_clade <- rowSums(presence[, in_clade, drop = FALSE])
count_other <- if (any(in_other)) rowSums(presence[, in_other, drop = FALSE]) else 0

# Keep drivers: >=3 in clade AND 0 elsewhere

keep_rows <- (count_clade >= 3) & (count_other == 0)

filtered_effect_mat <- effect_matrix_after_check_T[keep_rows, , drop = FALSE]


# limit columns to clade samples (ordered as given 

cols_to_show <- intersect(clade_samples, colnames(filtered_effect_mat))

# Rebuild top annotation for the subset of columns

ordered_subset <- cols_to_show
sample_group_sub <- ifelse(ordered_subset %in% P3G6_samples_filtered, "Ascites", NA)
cell_type_sub    <- unname(sample_to_cat[ordered_subset])

top_annot_sub <- HeatmapAnnotation(
  Location = sample_group_sub,
  Sample   = cell_type_sub,
  col = list(
    Location = group_colors,
    Sample   = celltype_colors
  ),
  annotation_name_side = "left"
)

# Rebuild row annotation aligned to filtered rows

driver_class_sub <- driver_class[match(rownames(filtered_effect_mat),
                                       rownames(effect_matrix_after_check_T))]
row_annot_sub <- rowAnnotation(
  DriverClass = driver_class_sub,
  col = list(DriverClass = driver_class_cols),
  annotation_name_side = "top"
)

oncoPrint(
  filtered_effect_mat[, ordered_subset, drop = FALSE],
  alter_fun              = alter_fun,
  col                    = col,
  column_title           = "P3G6 oncoplot (clade-only drivers ≥2)",
  heatmap_legend_param   = heatmap_legend_param,
  show_column_names      = TRUE,
  column_order           = ordered_subset,
  top_annotation         = top_annot_sub,
  left_annotation        = row_annot_sub,
  right_annotation       = NULL,
  remove_empty_rows      = TRUE,
  remove_empty_columns   = FALSE,
  column_names_gp        = gpar(fontsize = 10),
  row_names_gp           = gpar(fontsize = 10)
)


