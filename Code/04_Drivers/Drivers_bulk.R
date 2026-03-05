################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Driver analysis bulk
# Author: Alexander Steemers
################################################################################

# load libraries
library(dplyr)
library(stringr)
library(tibble)
library(stringi)
library(readxl)
library(VariantAnnotation)
library(GenomicRanges)
library(S4Vectors)
library(tidyr)
library(maftools)
library(writexl)

# load csv files that contain all somatic mutations
# these can be found in Mendeley doi: 10.17632/phk5jhhm7d.1 --> Drivers_bulk
# there should be 21 samples
# make a list called: all_csv_files

# load metadata
bulk_wgs_meta <- read_xlsx(Data/Bulk_sample_manuscript.xlsx)

# Read + combine
mutation_data_list <- lapply(all_csv_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  name <- tools::file_path_sans_ext(basename(f))
  barcode <- sub("_mutation_summary$", "", name)
  df$Tumor_Sample_Barcode <- barcode
  df
})

# generate one large dataframe
mutation_data <- bind_rows(mutation_data_list)

# trim df; only need first 38 columns (i.e. up to and including am_pathogenicity)
mutation_data <- mutation_data[, 1:38]

# Now start downstream filtering

# Strict filtering pipeline: 
# 1) keep only PASS
# 2) keep only HIGH/MODERATE
# 3) for missense mutations keep only likely_pathogenic
# 4) keep only protein coding variants
# 5) keep only Burkitt genes
# 6) remove duplicates

# 1) keep only PASS
mutation_data_PASS <- mutation_data %>% filter(FILTER == "PASS")

# 2) keep only HIGH/MODERATE
mutation_data_high_mod <- mutation_data_PASS %>%
  filter(IMPACT %in% c("HIGH", "MODERATE"))

mutation_data_high_mod$Consequence <- sub("&.*", "", mutation_data_high_mod$Consequence) # Make consequence annotations more simple

unique(mutation_data_high_mod$Consequence)

# Map to MAF-compatible classifications
map_to_maf_class <- function(consequence) {
  consequence <- tolower(consequence)
  
  if (grepl("frameshift_variant", consequence)) return("Frame_Shift")
  if (grepl("inframe_deletion|inframe_insertion", consequence)) return("Inframe_INDEL")
  if (grepl("stop_gained", consequence)) return("Nonsense_Mutation")
  if (grepl("stop_lost", consequence)) return("Nonstop_Mutation")
  if (grepl("missense_variant", consequence)) return("Missense_Mutation")
  # All else becomes "Other"
  return("Other")
}

mutation_data_high_mod$Variant_Classification <- sapply(mutation_data_high_mod$Consequence, map_to_maf_class)

# 3) for missense mutations keep only likely_pathogenic
mutation_data_filt <- mutation_data_high_mod %>% 
  filter( Consequence != "missense_variant" | (!am_class %in% c("likely_benign")))

# 4) keep only protein coding variants
mutation_data_filt_2 <- mutation_data_filt %>%
  filter(BIOTYPE == "protein_coding")

# Remove all non-coding mutations + synonymous
unique(mutation_data_filt_2$Consequence)

# Get known BL driver genes
BL_drivers <- read.delim(
  "Data/BL_Intogen_Panea_Blood2019_Nicole_Blood_2023_modified.txt",
  header = FALSE, stringsAsFactors = FALSE
)[,1]

# 5) keep only Burkitt genes
mutation_data_filt_3 <- mutation_data_filt_2[mutation_data_filt_2$SYMBOL %in% BL_drivers, ] # Filter based on known driver list

mutation_data_filt_3$Variant_Type <- ifelse(
  nchar(mutation_data_filt_3$REF) == 1 & nchar(mutation_data_filt_3$ALT) == 1,
  "SNP", "INDEL"
)

# 7) remove duplicates
mutation_data_filt_4 <- as.data.frame(mutation_data_filt_3) %>% 
  as_tibble() %>%
  mutate(has_am = !is.na(am_class) & nzchar(am_class)) %>%
  group_by(across(1:18)) %>%
  arrange(desc(has_am), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(-has_am)

mutation_data_filt_5 <- mutation_data_filt_4 %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(
    impact_chr = as.character(IMPACT),
    prefer_high = !is.na(impact_chr) & impact_chr == "HIGH"
  ) %>%
  group_by(across(1:14)) %>%
  arrange(desc(prefer_high), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  dplyr::select(-prefer_high, -impact_chr)

# manually check missense mutations that do not have alphamissense prediction in alphamissense database: https://alphamissense.hegelab.org/results
manual_check <- read_excel("Data/bulk_alphamissense_manual_check.xlsx") %>%
  mutate(
    Protein_position = as.character(Protein_position),
    Codons = as.character(Codons)
  )


keys <- c("SYMBOL", "Feature", "Protein_position", "Amino_acids",
          "Codons", "Variant_Classification", "Variant_Type")

# Make sure types are consistent (Excel often breaks these)
mutation_data_filt_5 <- mutation_data_filt_5 %>%
  mutate(
    Protein_position = as.character(Protein_position),
    Codons = as.character(Codons),
    am_pathogenicity = readr::parse_number(as.character(am_pathogenicity))
  )

manual_lookup <- manual_check %>%
  mutate(
    Protein_position = as.character(Protein_position),
    Codons = as.character(Codons),
    am_pathogenicity = readr::parse_number(as.character(am_pathogenicity))
  ) %>%
  group_by(across(all_of(keys))) %>%
  summarise(
    am_class = if (all(is.na(am_class))) NA_character_ else first(na.omit(am_class)),
    am_pathogenicity = if (all(is.na(am_pathogenicity))) NA_real_ else first(na.omit(am_pathogenicity)),
    .groups = "drop"
  )

# Build key strings for matching (no join, no row expansion)
key_x <- do.call(paste, c(mutation_data_filt_5[keys], sep = "||"))
key_y <- do.call(paste, c(manual_lookup[keys], sep = "||"))

idx <- match(key_x, key_y)   # positions in manual_lookup for each row in mutation_data_filt_5

# Fill only where we have a match (and optionally only where missing)
hit <- !is.na(idx)

mutation_data_filt_5$am_class[hit] <-
  ifelse(is.na(mutation_data_filt_5$am_class[hit]),
         manual_lookup$am_class[idx[hit]],
         mutation_data_filt_5$am_class[hit])

mutation_data_filt_5$am_pathogenicity[hit] <-
  ifelse(is.na(mutation_data_filt_5$am_pathogenicity[hit]),
         manual_lookup$am_pathogenicity[idx[hit]],
         mutation_data_filt_5$am_pathogenicity[hit])

mutation_data_filt_6 <- mutation_data_filt_5 %>%
  filter(
    !str_detect(Consequence, "missense_variant") |
      am_class %in% c("likely_pathogenic", "ambiguous")
  )

# save the data# save the data# save the dataframe to use in single-cell driver analysis
write.csv(mutation_data_filt_6, "mutation_data_filtered.csv", row.names = FALSE)

# Build maf_data with required columns
maf_data <- mutation_data_filt_6 %>%
  dplyr::select(
    Hugo_Symbol = SYMBOL,
    Chromosome = CHROM,
    Start_Position = POS,
    End_Position = POS,
    Reference_Allele = REF,
    Tumor_Seq_Allele2 = ALT,
    Tumor_Sample_Barcode = Sample,
    Variant_Classification,
    Variant_Type
  )

# Write to MAF file
maf_temp_file <- "bulk_mutations.maf"
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
  Sex = c("Female" = "#1D3557", "Male" = "#A9B7C8"),
  scWGS = c("Yes" = "#D2BD96", "No" = "#F0E9DC"),
  Translocation = c("MYC-IGH" = "#6C5B7B", "MYC-IGK" = "#C6BFD0")
)

# Generate oncoplot
pdf("oncoplot_bulk_rebuttal.pdf", width = 5.3, height = 10) 
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

genes_per_sample <- maf_object@data %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(num_genes = n_distinct(Hugo_Symbol))

# View summary statistics
summary(genes_per_sample$num_genes)

write_xlsx(maf_data, path = "maf_data_bulk_v3.xlsx")

