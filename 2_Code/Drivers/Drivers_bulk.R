################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look drivers in bulk tumour samples and generate oncoplot
# Author: Alexander Steemers
# Date: July 2025
################################################################################

# Load required libraries
library(VariantAnnotation)
library(dplyr)
library(maftools)
library(readxl)

# Set output and input base directories
setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Drivers/")
base_dir <- "~/hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus"
output_dir <- "~/hpc/pmc_vanboxtel/projects/Burkitt/2_Code/Drivers/Markus/Bulk_drivers"
base_dir   <- path.expand(base_dir)
output_dir <- path.expand(output_dir)

# Import metadata (optional)
bulk_wgs_meta <- read_xlsx(
  "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_manuscript.xlsx"
)

# Explicit list of VCF files
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

# --- CSQ helpers (robust to VEP field order and NA) ----
get_csq_indices <- function(vcf) {
  info_h <- VariantAnnotation::info(header(vcf))
  if (!"CSQ" %in% rownames(info_h))
    return(list(symbol_idx = NA_integer_, cons_idx = NA_integer_))
  desc <- info_h["CSQ", "Description"]
  fmt  <- sub(".*[Ff]ormat:\\s*", "", desc)
  fields <- strsplit(fmt, "\\|")[[1]]
  symbol_idx <- which(fields %in% c("SYMBOL","Gene"))[1]
  cons_idx   <- which(fields %in% c("Consequence"))[1]
  list(symbol_idx = ifelse(length(symbol_idx)==1, symbol_idx, NA_integer_),
       cons_idx   = ifelse(length(cons_idx)==1,   cons_idx,   NA_integer_))
}

extract_first_field <- function(csq_entry, target_idx) {
  if (is.na(target_idx)) return(NA_character_)
  if (length(csq_entry) == 0) return(NA_character_)
  if (length(csq_entry) > 1) csq_entry <- paste0(csq_entry, collapse = ",")
  if (is.na(csq_entry) || identical(csq_entry, "")) return(NA_character_)
  entries <- strsplit(csq_entry, ",", fixed = TRUE)[[1]]
  for (entry in entries) {
    fields <- strsplit(entry, "|", fixed = TRUE)[[1]]
    if (length(fields) >= target_idx) {
      val <- fields[[target_idx]]
      if (!is.null(val) && !is.na(val) && nzchar(val)) return(val)
    }
  }
  NA_character_
}

summarize_vcf <- function(vcf_path, out_dir = output_dir) {
  vcf <- VariantAnnotation::readVcf(vcf_path)
  rr  <- rowRanges(vcf)
  csq <- VariantAnnotation::info(vcf)$CSQ
  idx <- get_csq_indices(vcf)
  
  gene_symbols <- vapply(as.list(csq), extract_first_field, character(1), target_idx = idx$symbol_idx)
  consequences <- vapply(as.list(csq), extract_first_field, character(1), target_idx = idx$cons_idx)
  
  REF_vec <- as.character(mcols(rr)$REF)
  ALT_vec <- vapply(mcols(rr)$ALT, function(x) paste(as.character(x), collapse = ","), character(1))
  
  sample_label <- sub("\\.vep.*$", "", basename(vcf_path))
  
  df <- data.frame(
    Sample       = sample_label,
    File         = basename(vcf_path),
    seqnames     = as.character(seqnames(rr)),
    position     = start(rr),
    End_Position = end(rr),
    REF          = REF_vec,
    ALT          = ALT_vec,
    QUAL         = mcols(rr)$QUAL,
    FILTER       = as.character(mcols(rr)$FILTER),
    GENE         = gene_symbols,
    CONSEQUENCE  = consequences,
    stringsAsFactors = FALSE
  )
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  label <- tools::file_path_sans_ext(basename(vcf_path))
  outfile <- file.path(out_dir, paste0(label, "_mutation_summary.csv"))
  write.csv(df, outfile, row.names = FALSE)
  message("  [OK] Wrote: ", outfile)
}

# ------------------------ Main ------------------------
for (vcf_file in vcf_files) {
  # sample_name = filename before first dot
  sample_name <- strsplit(basename(vcf_file), "\\.", fixed = FALSE)[[1]][1]
  # folder = sample_name with everything after first underscore stripped
  # e.g. P3G6_bulk -> P3G6 ; PMCID211AAO_2 -> PMCID211AAO ; P856_D -> P856
  folder <- sub("_.*$", "", sample_name)
  
  vcf_path <- file.path(base_dir, folder, vcf_file)
  
  if (!file.exists(vcf_path)) {
    message("File not found: ", vcf_path, " — skipping.")
    next
  }
  
  message("Processing: ", vcf_path)
  tryCatch(
    summarize_vcf(vcf_path, out_dir = output_dir),
    error = function(e) message("  [WARN] Failed on ", vcf_file, ": ", conditionMessage(e))
  )
}

# Collect all mutation summary CSVs (written to output_dir)
all_csv_files <- list.files(
  path = output_dir,
  pattern = "_mutation_summary\\.csv$",
  full.names = TRUE,
  recursive = FALSE
)
message("Found ", length(all_csv_files), " mutation summary CSVs in: ", output_dir)

all_csv_files <- all_csv_files[!grepl("PVA9|PMCID690AAT", all_csv_files)] # the first two are removed because I already have bulk of those and the third is removed because no MYC translocation was found in that sample

# Read each CSV and assign a unique Tumor_Sample_Barcode based on filename
mutation_data_list <- lapply(all_csv_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  # Extract something like "P856_BCELLBULK" from the filename
  name <- tools::file_path_sans_ext(basename(f))
  barcode <- sub("_mutation_summary$", "", name)
  
  df$Tumor_Sample_Barcode <- barcode
  df
})

# Combine everything into one data frame
mutation_data <- do.call(rbind, mutation_data_list)

# Simplify CONSEQUENCE column by taking only the part before "&"

mutation_data$CONSEQUENCE_ <- sub("&.*", "", mutation_data$CONSEQUENCE)

# Remove all non-coding mutations + synonymous

unique(mutation_data$CONSEQUENCE)

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

mutation_data_filt <- mutation_data[grepl(pattern, mutation_data$CONSEQUENCE, ignore.case = TRUE), ]

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

mutation_data_filt_2 <- mutation_data_filt[mutation_data_filt$GENE %in% BL_drivers, ] # Filter based on known driver list

mutation_data_filt_3 <- mutation_data_filt_2[mutation_data_filt_2$FILTER != "SnpCluster", ] # Filter out SnpClusters

mutation_data_filt_3$CONSEQUENCE <- sub("&.*", "", mutation_data_filt_3$CONSEQUENCE) # Make consequence annotations more simple

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
mutation_data_filt_3$Variant_Classification <- sapply(mutation_data_filt_3$CONSEQUENCE, map_to_maf_class)

mutation_data_filt_3$Variant_Type <- ifelse(
  nchar(mutation_data_filt_3$REF) == 1 & nchar(mutation_data_filt_3$ALT) == 1,
  "SNP", "INDEL"
)

# Build maf_data with required columns
maf_data <- mutation_data_filt_3 %>%
  dplyr::select(
    Hugo_Symbol = GENE,
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
pdf("Figures/oncoplot_bulk_wo_690AAT.pdf", width = 5.3, height = 10) 
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

