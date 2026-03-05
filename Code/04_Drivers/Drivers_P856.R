################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Script to plot oncoplot for P856
# Author: Alexander Steemers
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
library(stringr)
library(writexl)

# Find blacklisted samples

low_callable_df<- read.csv('low_callable_loci.csv') # generated after running 03_QC_step1.R
below_curve_df <- read.csv('below_curve_samples.csv') # generated after running 03_QC_step1.R
bad_baf_df     <- read.csv('bad_baf_samples.csv') # generated after running 03_QC_step1.R
fail_vaf_df    <- read.csv('PTA_samples_failVAFcheck.txt') # generated after running 04_QC_step2.R

blacklist <- unique(c(below_curve_df$Sample_name,
                      low_callable_df$Sample_name,
                      bad_baf_df$Sample_name,
                      fail_vaf_df$samplename))

# Desired column order (from CellPhy tree)
P856_samples <-  c( "P856GDDUBC42", "P856GDDUBC44", "PB14458-BLPL-BCELLP4L3", "PB14458-BLPL-BCELLP4K5", "PB14458-BLPL-BCELLP4M3", "PB14458-BLPL-BCELLP4B3", "PB14458-BLPL-BCELLP4L5" , "P856GDDBBC63", "PB14458-BLBM-BCELLP2F2", "PB14458-BLBM-BCELLP2B3", "P856GDDBBC60", "P856GDDBBC46", "PB14458-BLBM-BCELLP2B4", "PB14458-BLBM-BCELLP2C4", "PB14458-BLBM-BCELLP2N2", "P856GDDBBC59", "P856GDDBBC64", "PB14458-BLBM-BCELLP2F4", "P856GDDBBC58", "P856GDDBBC54", "P856GDDBBC48", "PB14458-BLBM-BCELLP2I2", "PB14458-BLBM-BCELLP2L4",  "PB14458-BLBM-BCELLP2E4", "P856GDDBBC57", "P856GDDBBC62" ,"PB14458-BLBM-BCELLP2L3", "P856GDDBBC61","PB14458-BLPL-BCELLBULK", "PB14458-BLBM-BCELLBULK") 
P856_normal <- c( "P856GDDUBC42", "P856GDDUBC44")

P856_tumour <- c("PB14458-BLPL-BCELLP4L3", "PB14458-BLPL-BCELLP4K5", "PB14458-BLPL-BCELLP4M3", "PB14458-BLPL-BCELLP4B3", "PB14458-BLPL-BCELLP4L5" , "P856GDDBBC63", "PB14458-BLBM-BCELLP2F2", "PB14458-BLBM-BCELLP2B3", "P856GDDBBC60", "P856GDDBBC46", "PB14458-BLBM-BCELLP2B4", "PB14458-BLBM-BCELLP2C4", "PB14458-BLBM-BCELLP2N2", "P856GDDBBC59", "P856GDDBBC64", "PB14458-BLBM-BCELLP2F4", "P856GDDBBC58", "P856GDDBBC54", "P856GDDBBC48", "PB14458-BLBM-BCELLP2I2", "PB14458-BLBM-BCELLP2L4",  "PB14458-BLBM-BCELLP2E4", "P856GDDBBC57", "P856GDDBBC62" ,"PB14458-BLBM-BCELLP2L3", "P856GDDBBC61","PB14458-BLPL-BCELLBULK", "PB14458-BLBM-BCELLBULK")

# Get known BL driver genes

BL_drivers <- read.delim(
  "Data/BL_Intogen_Panea_Blood2019_Nicole_Blood_2023_modified.txt",
  header = FALSE, stringsAsFactors = FALSE
)[,1]

# Function to extract driver mutations

extract_genes_mutations <- function(vcf) {
  # --- 1) Parse CSQ field layout from header ---
  csq_field  <- VariantAnnotation::info(vcf)$CSQ
  csq_desc   <- as.character(VariantAnnotation::info(header(vcf))["CSQ", "Description"])
  csq_format <- sub(".*Format:\\s*", "", csq_desc)
  csq_fields <- strsplit(csq_format, "\\|")[[1]]
  
  # indices of useful fields
  gene_idx   <- which(csq_fields %in% c("SYMBOL","Gene","Gene_Name"))[1]
  conseq_idx <- which(csq_fields == "Consequence")[1]
  hgvsp_idx  <- which(csq_fields == "HGVSp")[1]
  hgvsc_idx  <- which(csq_fields == "HGVSc")[1]
  am_class_idx <- which(csq_fields == "am_class")[1]
  am_path_idx  <- which(csq_fields == "am_pathogenicity")[1]
  protpos_idx <- which(csq_fields %in% c("Protein_position","ProteinPosition"))[1]
  aa_idx      <- which(csq_fields %in% c("Amino_acids","AminoAcids"))[1]
  
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
  parsed     <- strsplit(csq_flat, "\\|", fixed = FALSE)
  max_len    <- max(vapply(parsed, length, 0L))
  parsed_mat <- t(vapply(parsed, function(x){ length(x) <- max_len; x }, character(max_len)))
  
  ann_df <- data.frame(
    variant     = variant_ids[var_index],
    chrom       = chrom[var_index],
    pos         = pos[var_index],
    ref         = REF[var_index],
    alt         = ALT[var_index],
    gene        = if (!is.na(gene_idx))   parsed_mat[, gene_idx]   else NA,
    consequence = if (!is.na(conseq_idx)) parsed_mat[, conseq_idx] else NA,
    protein_position = if (!is.na(protpos_idx)) parsed_mat[, protpos_idx] else NA,
    amino_acids      = if (!is.na(aa_idx))      parsed_mat[, aa_idx]      else NA,
    hgvsp       = if (!is.na(hgvsp_idx))  parsed_mat[, hgvsp_idx]  else NA,
    hgvsc       = if (!is.na(hgvsc_idx))  parsed_mat[, hgvsc_idx]  else NA,
    am_class         = if (!is.na(am_class_idx)) parsed_mat[, am_class_idx] else NA,
    am_pathogenicity = if (!is.na(am_path_idx))  parsed_mat[, am_path_idx]  else NA,
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
                         "protein_position","amino_acids",
                         "am_class","am_pathogenicity",
                         "sample","ref_count","alt_count","dp","vaf") ])
  out
}

# obtain vcfs from Mendeley doi: 10.17632/phk5jhhm7d.1 --> Drivers_single_cell --> P856
# make list called vcf_files

# read all VCFs, extract tables, and combine

genes_mut_merged <- vcf_files %>%
  map(~ VariantAnnotation::readVcf(.x)) %>%
  map(extract_genes_mutations) %>%
  bind_rows()

# Remove duplicate rows

genes_mut_merged_nodup <- unique(genes_mut_merged)

# remove low vaf
genes_mut_merged_highvaf <- genes_mut_merged_nodup %>%
  filter(vaf >= 0.15)

# remove normal B-cells
genes_mut_merged_filterd_1 <- genes_mut_merged_highvaf %>%
  filter(!sample %in% P856_normal)

# Remove rows with samples which have been blacklisted
genes_mut_merged_filterd_2 <- genes_mut_merged_filterd_1 %>%
  filter(!sample %in% blacklist)

# Filter out rows which do not have genes in the BL list of known genes

genes_mut_merged_filterd_3 <- genes_mut_merged_filterd_2 %>% filter(gene %in% BL_drivers)

# Simplify mutation annotation

genes_mut_merged_filterd_4 <- genes_mut_merged_filterd_3 %>%
  mutate(consequence = stringr::str_squish(stringr::str_replace(consequence, "&.*$", "")))

# keep only likely_pathogenic for misssense mutations
genes_mut_merged_filterd_5 <- genes_mut_merged_filterd_4 %>%
  filter(
    consequence != "missense_variant" |
      am_class %in% c("likely_pathogenic", "ambiguous")
  )

# Remove rows with VAF 0 or N/A

genes_mut_merged_filterd_6 <- genes_mut_merged_filterd_5 %>%
  filter(!is.na(vaf), vaf != 0)

# Keep only protein coding variants

driver_terms <- c(
  "stop_gained", "inframe_insertion", "missense_variant", "frameshift_variant", "splice_donor_variant"    
  ,"inframe_deletion","protein_altering_variant", "splice_acceptor_variant", "start_lost", "stop_lost"     
)

genes_mut_merged_filterd_7 <- subset(genes_mut_merged_filterd_6, grepl(paste(driver_terms, collapse="|"), consequence))

# Provide sample info

input_df <-  read_excel('Data/Sample_overview.xlsx') 
P856_rows <- input_df %>% filter(Novogene_ID == "P856") %>% arrange(factor(Myc_translocation_IGV, levels = c("no", "yes")))   
normal_bcell_samples <- P856_rows$Sample_name[P856_rows$Myc_translocation_IGV == "No"]

# Remove duplicates
genes_mut_merged_filterd_7dedup <- genes_mut_merged_filterd_7 %>%
  distinct()


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

genes_mut_merged_filterd_11 <- genes_mut_merged_filterd_11 %>%
  mutate(am_pathogenicity = as.numeric(am_pathogenicity))

# look at bulk sample
bulk_samples <- read.csv("mutation_data_filtered.csv") # generated after running Drivers_bulk.R script
bulk_sample <- bulk_samples[bulk_samples$Sample %in% c("P856_D", "P856_R"), ]
bulk_sample <- bulk_sample %>%
  dplyr::rename(sample = Sample) %>%                 # rename column
  mutate(sample = ifelse(sample == "P856_D", "PB14458-BLPL-BCELLBULK", sample)) %>% 
  mutate(sample = ifelse(sample == "P856_R", "PB14458-BLBM-BCELLBULK", sample)) %>% 
  mutate(driver = paste0(SYMBOL, "_", CHROM, ":", POS, "_", REF, ">", ALT)) %>%
  dplyr::rename(consequence = Consequence)

combined_data <- bind_rows(
  genes_mut_merged_filterd_11,
  bulk_sample
)

# save combined_data 
write.csv(combined_data, "combined_data_P856.csv")

# Prepare oncoplot matrix 

oncoplot_matrix <- combined_data %>%
  dplyr::select(driver, sample, consequence) %>%      # keep only relevant columns
  pivot_wider(
    names_from = sample,
    values_from = consequence,
    values_fill = list(mutation = NA)        # fill missing with NA
  )


# Ensure all P856_samples exist as columns (add NA if missing)

for (s in P856_samples) {
  if (!s %in% colnames(oncoplot_matrix)) {
    oncoplot_matrix[[s]] <- NA
  }
}

oncoplot_matrix <- oncoplot_matrix %>%
  dplyr::select(-any_of(P856_normal))

# Reorder columns: driver first, then P856_samples
oncoplot_matrix <- oncoplot_matrix %>%
  dplyr::select(driver, all_of(P856_tumour))

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
P856_samples_filtered <- setdiff(P856_samples, blacklist)
P856_samples_filtered <- intersect(P856_samples_filtered, colnames(effect_matrix_after_check_T))

# Cell type colors and vector aligned to ordered_samples
celltype_colors <- c("Normal B-cell" = "#E7872B",
                     "Burkitt Lymphoma cell" = "#3F78C1",
                     "Bulk Tumour" = "#00008B")

# Categorize samples in P856_rows
P856_rows <- P856_rows %>%
  mutate(
    Category = case_when(
      Sample_name %in% normal_bcell_samples ~ "Normal B-cell",
      Sample_name == "PB14458-BLPL-BCELLBULK" ~ "Bulk Tumour",
      Sample_name == "PB14458-BLBM-BCELLBULK" ~ "Bulk Tumour",
      TRUE ~ "Burkitt Lymphoma cell"
    ),
    Category = factor(Category,
                      levels = c("Normal B-cell", "Burkitt Lymphoma cell", "Bulk Tumour"))
  )

# Make a named lookup and align it to ordered_samples
sample_to_cat <- setNames(as.character(P856_rows$Category), P856_rows$Sample_name)
cell_type <- unname(sample_to_cat[ordered_samples])

# Variant class colors (ensure these strings match your 'mutation' values)
col <- c(
  missense_variant    = "#008000",   # green
  intron_variant      = "#7D3C98",   # purple
  inframe_insertion   = "#8B0000",   # dark red
  frameshift_variant  = "#1F77B4",   # blue
  protein_altering_variant         = "#FF8C00",  
  stop_gained         = "#FF0000",   # red
  inframe_deletion    = "#8B4513",   # brown
  start_lost = "#A9A9A9",   # dark gray
  splice_donor_variant            = "#40E0D0",   # turquoise
  splice_acceptor_variant              = "#00CED1"   # dark turquoise
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

column_title <- "P856 oncoplot"

# Plot
pdf("P856_oncoprint_rebuttal.pdf",
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
