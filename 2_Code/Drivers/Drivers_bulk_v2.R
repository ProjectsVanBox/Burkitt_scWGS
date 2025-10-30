################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Look for driver alterations in bulk tumour samples and generate oncoplot
# Author: Alexander Steemers
# Date: July 2025
################################################################################

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(S4Vectors)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(maftools)
  library(writexl)
})

# --------------------------- Configuration ------------------------------------

# Set working/output dirs
setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/")
base_dir  <- "~/hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus"
output_dir <- "~/hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus/Bulk_drivers"

base_dir   <- path.expand(base_dir)
output_dir <- path.expand(output_dir)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Optional metadata (not strictly required for oncoplot)
meta_path <- "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_manuscript.xlsx"
if (file.exists(meta_path)) {
  bulk_wgs_meta <- read_xlsx(meta_path)
}

# Explicit list of VCFs (relative to <base_dir>/<folder>/file)
vcf_files <- c(
  "PMCID104AAO.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID132AAL.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID137AAO.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID163AAN.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID211AAO_2.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID321AAO.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID340AAO.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID458AAQ.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID491AAS.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID509AAT.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID540AAN.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID610AAS.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID690AAT.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID821AAL.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID867AAT.vep.effect.somatic.pass_snpcluster.vcf",
  "PMCID967AAP.vep.effect.somatic.pass_snpcluster.vcf",
  "P3G6_bulk.vep.effect.somatic.pass_snpcluster.vcf",
  "PRN4_bulk.vep.effect.somatic.pass_snpcluster.vcf",
  "PIA9_bulk.vep.effect.somatic.pass_snpcluster.vcf",
  "PVA9_bulk.vep.effect.somatic.pass_snpcluster.vcf",
  "PJBU_bulk.vep.effect.somatic.pass_snpcluster.vcf",
  "P856_D_vep.effect.somatic.pass_snpcluster.vcf",
  "P856_R_vep.effect.somatic.pass_snpcluster.vcf"
)

# Exclusions (optional; you had these earlier)
exclude_patterns <- "PVA9|PMCID690AAT"

# -------------------------- VEP CSQ parsing utils -----------------------------

`%||%` <- function(a, b) if (!is.na(a)) a else b

get_csq_map <- function(vcf) {
  info_h <- VariantAnnotation::info(header(vcf))
  if (!"CSQ" %in% rownames(info_h)) {
    stop("No CSQ field found in VCF header.")
  }
  desc <- info_h["CSQ", "Description"]
  fmt  <- sub(".*[Ff]ormat:\\s*", "", desc)
  fields <- strsplit(fmt, "\\|")[[1]]
  idx <- function(name) { w <- which(fields == name); if (length(w)) w[1] else NA_integer_ }
  list(
    fields       = fields,
    i.Allele      = idx("Allele"),
    i.Consequence = idx("Consequence"),
    i.IMPACT      = idx("IMPACT"),
    i.SYMBOL      = idx("SYMBOL"),
    i.Gene        = idx("Gene"),
    i.BIOTYPE     = idx("BIOTYPE")
  )
}

# Parse CSQ value (which can be length 0/1/>1) into a tibble of rows (one per CSQ entry)
parse_csq_entries <- function(csq_value, map) {
  if (length(csq_value) == 0 || all(is.na(csq_value))) return(NULL)
  csq_str <- paste0(csq_value[!is.na(csq_value) & nzchar(csq_value)], collapse = ",")
  if (!nzchar(csq_str)) return(NULL)
  
  entries <- strsplit(csq_str, ",", fixed = TRUE)[[1]]
  out <- lapply(entries, function(entry) {
    f <- strsplit(entry, "|", fixed = TRUE)[[1]]
    get <- function(i) if (!is.na(i) && length(f) >= i && nzchar(f[i])) f[i] else NA_character_
    tibble::tibble(
      SYMBOL       = get(map$i.SYMBOL),
      Gene         = get(map$i.Gene),
      Consequence  = get(map$i.Consequence),
      IMPACT       = get(map$i.IMPACT),
      BIOTYPE      = get(map$i.BIOTYPE)
    )
  })
  dplyr::bind_rows(out)
}

# Summarize a VCF into long format: one row per (variant x gene consequence)
summarize_vcf_long <- function(vcf_path) {
  vcf <- VariantAnnotation::readVcf(vcf_path)
  rr  <- rowRanges(vcf)
  info_list <- VariantAnnotation::info(vcf)
  
  if (!"CSQ" %in% colnames(info_list)) {
    message("[WARN] No CSQ in: ", basename(vcf_path))
    return(tibble::tibble())
  }
  csq <- info_list$CSQ
  map <- get_csq_map(vcf)
  
  REF_vec <- as.character(mcols(rr)$REF)
  ALT_vec_full <- mcols(rr)$ALT
  ALT_vec <- vapply(ALT_vec_full, function(x) paste(as.character(x), collapse=","), character(1))
  
  sample_label <- sub("\\.vep.*$|_vep.*$|\\.vcf(?:\\.gz)?$", "", basename(vcf_path))
  
  rows <- vector("list", length(csq))
  for (i in seq_along(csq)) {
    ann <- parse_csq_entries(csq[[i]], map)
    if (is.null(ann)) next
    rows[[i]] <- ann %>%
      mutate(
        Sample       = sample_label,
        File         = basename(vcf_path),
        seqnames     = as.character(seqnames(rr))[i],
        position     = start(rr)[i],
        End_Position = end(rr)[i],
        REF          = REF_vec[i],
        ALT          = ALT_vec[i],
        QUAL         = mcols(rr)$QUAL[i],
        FILTER       = as.character(mcols(rr)$FILTER[i])
      )
  }
  
  df <- bind_rows(rows) %>%
    filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    relocate(Sample, File, seqnames, position, End_Position, REF, ALT, QUAL, FILTER,
             SYMBOL, Consequence, IMPACT, BIOTYPE, .before = Gene)
  
  df
}

# Writer wrapper – save one CSV per input VCF
summarize_vcf_write <- function(vcf_path, out_dir) {
  df <- summarize_vcf_long(vcf_path)
  if (nrow(df) == 0L) {
    message("  [WARN] No usable CSQ rows for: ", basename(vcf_path))
    return(invisible(NULL))
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  label <- tools::file_path_sans_ext(basename(vcf_path))
  outfile <- file.path(out_dir, paste0(label, "_mutation_summary.csv"))
  write.csv(df, outfile, row.names = FALSE)
  message("  [OK] Wrote: ", outfile)
  invisible(outfile)
}

# --------------------------- Main processing ----------------------------------

# Process each VCF -> CSV summary
for (vcf_file in vcf_files) {
  # e.g. P3G6_bulk.vep... -> folder "P3G6"; PMCID211AAO_2 -> "PMCID211AAO"; P856_D -> "P856"
  sample_name <- strsplit(basename(vcf_file), "\\.", fixed = FALSE)[[1]][1]
  folder <- sub("_.*$", "", sample_name)
  vcf_path <- file.path(base_dir, folder, vcf_file)
  
  if (!file.exists(vcf_path)) {
    message("File not found: ", vcf_path, " — skipping.")
    next
  }
  message("Processing: ", vcf_path)
  tryCatch(
    summarize_vcf_write(vcf_path, out_dir = output_dir),
    error = function(e) message("  [WARN] Failed on ", vcf_file, ": ", conditionMessage(e))
  )
}

# Collect CSVs
all_csv_files <- list.files(output_dir, pattern = "_mutation_summary\\.csv$", full.names = TRUE)
message("Found ", length(all_csv_files), " mutation summary CSVs in: ", output_dir)

# Optional: drop some
if (nzchar(exclude_patterns)) {
  keep <- !grepl(exclude_patterns, all_csv_files)
  all_csv_files <- all_csv_files[keep]
  message("Keeping ", length(all_csv_files), " after exclusions: ", exclude_patterns)
}

# Read + combine
mutation_data_list <- lapply(all_csv_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  name <- tools::file_path_sans_ext(basename(f))
  barcode <- sub("_mutation_summary$", "", name)
  df$Tumor_Sample_Barcode <- barcode
  df
})

mutation_data <- bind_rows(mutation_data_list)

mutation_data_filt <- mutation_data %>%
  filter(BIOTYPE == "protein_coding")

# Simplify CONSEQUENCE column by taking only the part before "&"

mutation_data_filt$Consequence <- sub("&.*", "", mutation_data_filt$Consequence)

# Remove all non-coding mutations + synonymous

unique(mutation_data_filt$Consequence)

# pattern of coding terms to KEEP

keep_cons <- c(
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

pattern <- paste(keep_cons, collapse = "|")

mutation_data_filt_1 <- mutation_data_filt[grepl(pattern, mutation_data_filt$Consequence, ignore.case = TRUE), ]

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

mutation_data_filt_2 <- mutation_data_filt_1[mutation_data_filt_1$SYMBOL %in% BL_drivers, ] # Filter based on known driver list

mutation_data_filt_3 <- mutation_data_filt_2[mutation_data_filt_2$FILTER != "SnpCluster", ] # Filter out SnpClusters

mutation_data_filt_3$Consequence <- sub("&.*", "", mutation_data_filt_3$Consequence) # Make consequence annotations more simple

# Map to MAF-compatible classifications
map_to_maf_class <- function(consequence) {
  consequence <- tolower(consequence)
  
  if (grepl("frameshift_variant", consequence)) return("Frame_Shift")
  if (grepl("inframe_deletion|inframe_insertion", consequence)) return("Frame_Shift")
  if (grepl("stop_gained", consequence)) return("Nonsense_Mutation")
  if (grepl("stop_lost", consequence)) return("Nonstop_Mutation")
  if (grepl("missense_variant", consequence)) return("Missense_Mutation")
  # All else becomes "Other"
  return("Other")
}

# Apply classification and determine variant type
mutation_data_filt_3$Variant_Classification <- sapply(mutation_data_filt_3$Consequence, map_to_maf_class)

mutation_data_filt_3$Variant_Type <- ifelse(
  nchar(mutation_data_filt_3$REF) == 1 & nchar(mutation_data_filt_3$ALT) == 1,
  "SNP", "INDEL"
)

mutation_data_filt_3_nodup <- mutation_data_filt_3 %>%
  distinct()                # drop exact duplicate rows

# Build maf_data with required columns
maf_data <- mutation_data_filt_3_nodup %>%
  dplyr::select(
    Hugo_Symbol = SYMBOL,
    Chromosome = seqnames,
    Start_Position = position,
    End_Position,
    Reference_Allele = REF,
    Tumor_Seq_Allele2 = ALT,
    Tumor_Sample_Barcode,
    Variant_Classification,
    Variant_Type
  )

maf_data$Tumor_Sample_Barcode <- gsub(
  "_bulk\\.vep\\.effect\\.somatic\\.pass_snpcluster$", "",
  maf_data$Tumor_Sample_Barcode
)

maf_data$Tumor_Sample_Barcode <- gsub(
  "\\.vep\\.effect\\.somatic\\.pass_snpcluster$", "",
  maf_data$Tumor_Sample_Barcode
)

maf_data$Tumor_Sample_Barcode <- gsub(
  "\\_vep\\.effect\\.somatic\\.pass_snpcluster$", "",
  maf_data$Tumor_Sample_Barcode
)

maf_data$Tumor_Sample_Barcode <- gsub(
  "^PMCID", "",
  maf_data$Tumor_Sample_Barcode
)

maf_data$Tumor_Sample_Barcode <- gsub(
  "_2$", "",
  maf_data$Tumor_Sample_Barcode
)

# Add new IDs

maf_data <- maf_data %>%
  left_join(bulk_wgs_meta[, c("Tumor_Sample_Barcode", "ID")], by = "Tumor_Sample_Barcode")

maf_data <- maf_data %>%
  dplyr::select(-Tumor_Sample_Barcode) %>%
  dplyr::rename(Tumor_Sample_Barcode = ID)

# Manually checked PRN4 and P856
# Saw that SMARCA4 was present in PRN4_D1 but in low vaf
# Also, HLA-A was present in PRN4_D2 (with vaf close to 0.5) but was missed because of control sample having low percentage of tumour contamination
# Add these manually 

# ---- HLA-A: copy PRN4_D1 -> PRN4_D2 ----
hla_row <- maf_data[maf_data$Hugo_Symbol == "HLA-A" & 
                      maf_data$Tumor_Sample_Barcode == "PRN4_D1", ]
hla_copy <- hla_row
hla_copy$Tumor_Sample_Barcode <- "PRN4_D2"
maf_data <- rbind(maf_data, hla_copy)

# ---- SMARCA4: copy PRN4_D2 -> PRN4_D1 ----
smarca4_row <- maf_data[maf_data$Hugo_Symbol == "SMARCA4" & 
                          maf_data$Tumor_Sample_Barcode == "PRN4_D2", ]
smarca4_copy <- smarca4_row
smarca4_copy$Tumor_Sample_Barcode <- "PRN4_D1"
maf_data <- rbind(maf_data, smarca4_copy)

# Write to MAF file

maf_temp_file <- "Data/bulk_mutations.maf"
write.table(maf_data, maf_temp_file, sep = "\t", quote = FALSE, row.names = FALSE)

all_vcs <- unique(maf_data$Variant_Classification)

maf_object <- read.maf(
  maf = maf_temp_file,
  vc_nonSyn = all_vcs
)

# Get a base palette of 7 distinct colors
custom_palette_6 <- c(
  "#1F78B4",  
  "#E31A1C",  
  "#33A02C",
  "#f0c571",  
  "#360f5a",  
  "#A6CEE3"   
)

# Assign names matching your variant classifications
names(custom_palette_6) <- c(
  "Frame_Shift",
  "Nonsense_Mutation",
  "Missense_Mutation",
  "Nonstop_Mutation",
  "Multi_Hit",
  "Other"
)

# Subset relevant clinical data
clinical_data <- bulk_wgs_meta[, c( "ID", "Timepoint", "Sex", "scWGS", "Translocation")]
clinical_data <- distinct(clinical_data)
clinical_data <- dplyr::rename(clinical_data, Tumor_Sample_Barcode = ID)

# Add clinical data to the maf_object
maf_object@clinical.data <- dplyr::left_join(maf_object@clinical.data, clinical_data, by = "Tumor_Sample_Barcode")

ann_colors <- list(
  Timepoint = c("Diagnosis" = "#A62639", "Relapse" = "#E4A3AD"),
  Sex = c("Male" = "#1D3557", "Female" = "#A9B7C8"),
  scWGS = c("No" = "#D2BD96", "Yes" = "#F0E9DC"),
  Translocation = c("MYC-IGH" = "#6C5B7B", "MYC-IGK" = "#C6BFD0")
)

# Mark BL vs Cancer-only in the gene name itself
genes_in_maf <- unique(maf_object@data$Hugo_Symbol)

black_genes <- intersect(genes_in_maf, BL_drivers)                         # BL (incl. overlap)
grey_genes  <- intersect(genes_in_maf, setdiff(Cancer_drivers, BL_drivers))# Cancer-only

# Generate oncoplot
pdf("Figures/oncoplot_bulk_wo_690AAT_v2.pdf", width = 5.3, height = 10) 
oncoplot(
  maf = maf_object,
  top = 1000,
  showTumorSampleBarcodes = TRUE,
  fontSize = 0.7,
  SampleNamefontSize = 1,
  legend_height = 15,
  colors = custom_palette_6,
  clinicalFeatures = c("Timepoint", "Sex", "scWGS", "Translocation"),
  annotationColor = ann_colors
)
dev.off()

print(summary(tmb(maf = maf_object)$total))

genes_per_sample <- maf_object@data %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(num_genes = n_distinct(Hugo_Symbol))

# View summary statistics
summary(genes_per_sample$num_genes)

write_xlsx(maf_data, path = "Data/maf_data_bulk.xlsx")


