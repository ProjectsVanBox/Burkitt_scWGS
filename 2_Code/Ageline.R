################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to plot ageline of single cell WGS samples
# Author: Alexander Steemers
# Date: June 2025
# Modified: 
################################################################################

# Load libraries

library(BSgenome)
library(VariantAnnotation)
library(RColorBrewer)
library(tibble)
library(grid)
library(readxl)
library(Biostrings)
library(GenomicRanges)
library(stringr)
library(nlme)
library(dplyr)
library(ggrepel)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# Load colour palettes

mycols_paired <- brewer.pal(12,"Paired")
mycols_dark2 <- brewer.pal(8, "Dark2")

# Set directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Ageline/")

# Load functions and plotting functions

source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Set date

date <- format(Sys.Date(), "%Y%m%d")

# Load metadata 

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
input_df_sc <- input_df[input_df$ResolveDNA_version %in% c("v1", "v2", "v2.0"), ]
input_df_sub <- input_df_sc[!is.na(input_df_sc$Callable_fraction) & !is.na(input_df_sc$Mean_coverage),]

diagnostic_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_overview.csv')

below_curve_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv')
fail_vaf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')
low_callable_loci <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv')

# Make blacklist

blacklist <- unique(c(below_curve_df$Sample_name, fail_vaf_df$samplename, low_callable_loci$Sample_name))

# Load single cell PASS autosomal variants

SBSs_autosomal_PASS <- readRDS(file = "../MutLoad/Data/autosomal_PASS_variants_VAF015.RDS")

# Filter out blacklist samples and bulk sample

input_df_sub <- input_df_sub[ !input_df_sub$Sample_name %in% blacklist, ]
SBSs_autosomal_PASS <- SBSs_autosomal_PASS[ !(names(SBSs_autosomal_PASS) %in% blacklist) ]

# Make dataframe and merge with input_df

mut_load_df <- data.frame(
  Sample_name        = names(SBSs_autosomal_PASS),   
  Number_of_mutations = lengths(SBSs_autosomal_PASS), 
  row.names           = NULL,
  stringsAsFactors    = FALSE
)

merged_df <- input_df_sub %>% left_join(mut_load_df, by = "Sample_name")

# Normalize mutations to callable loci

hg38_autosomal_nonN_genome_size <- 2745186691
merged_df$Callable_Loci <- as.numeric(merged_df$Callable_Loci)

merged_df$Nmut_adj <- merged_df$Number_of_mutations/merged_df$Callable_Loci * hg38_autosomal_nonN_genome_size

# Name samples based on myc translocation and SBS9 

merged_df <- merged_df %>% mutate(MYC_SBS9_status = paste0(Myc_translocation_IGV, "_", SBS9))

merged_df <- merged_df %>% mutate(
    Sample = case_when(
      MYC_SBS9_status == "Yes_Positive" ~ "Burkitt Lymphoma cell",
      MYC_SBS9_status == "No_Positive"  ~ "Normal Memory B cell",
      MYC_SBS9_status == "No_Negative"  ~ "Normal Naive B cell",
      MYC_SBS9_status == "No_NA"  ~ "Unknown",
      MYC_SBS9_status == "Yes_NA"  ~ "Burkitt Lymphoma cell")
    )

# Get slopes and co-efficients

eq_df <- merged_df %>%
  group_by(Sample) %>%
  summarise(m = coef(lm(Nmut_adj ~ Age_at_sampling_Y))[2],
            c = coef(lm(Nmut_adj ~ Age_at_sampling_Y))[1],
            .groups = "drop") %>%
  mutate(label = paste0("y = ", round(m, 0), "x + ", round(c, 0)),
         x_pos = 20,                      
         y_pos = m * x_pos + c + 100)   

# Make ageline plot 

pdf("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Ageline/Figures/Ageline_lm.pdf",
    width = 10, height = 6)

ggplot(merged_df, aes(Age_at_sampling_Y, Nmut_adj, colour = Sample)) +
  #geom_smooth(data = subset(merged_df, Sample == "Burkitt Lymphoma cell"),
  #            method = "lm", se = FALSE, linewidth = 0.7,
  #            fullrange = TRUE, colour = "#3F78C1") +
  #geom_smooth(data = subset(merged_df, Sample == "Normal Memory B cell"),
  #            method = "lm", se = FALSE, linewidth = 0.7,
  #            fullrange = TRUE, colour = "#B96C22") +
  #geom_smooth(data = subset(merged_df, Sample == "Normal Naive B cell"),
  #            method = "lm", se = FALSE, linewidth = 0.7,
  #            fullrange = TRUE, colour = "#EC9F55") +
  #geom_smooth(data = subset(merged_df, Sample == "Unknown"),
  #            method = "lm", se = FALSE, linewidth = 0.7,
  #            fullrange = TRUE, colour = "#000000") +
  geom_jitter(width = 0.4, height = 0.4, size = 3) +
 
  #geom_text(data = eq_df,
  #          aes(x = x_pos, y = y_pos, label = label, colour = Sample),
  #          hjust = -0.05, vjust = 1, size = 5, show.legend = FALSE) +
  
  scale_colour_manual(values = c("Burkitt Lymphoma cell" = "#3F78C1",
                                 "Normal Memory B cell"  = "#B96C22",
                                 "Normal Naive B cell"   = "#EC9F55",
                                 "Unknown" = "#000000" )) +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 3500), breaks = seq(0,3500,500), expand = c(0, 0)) +
  labs(x = "Age at sampling (in years)", y = "SNVs per cell") +
  
  coord_cartesian(clip = "off") +   
  geom_abline(intercept = 1139, slope = 17, linetype = "dashed", colour = "black") +
  geom_abline(intercept = 215, slope = 15, linetype = "dotted", colour = "black") +
  theme_CHemALL() +
  theme(
    legend.position      = c(1.05, 1),   
    legend.justification = c(0, 1),      
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.margin        = margin(2, 4, 2, 4),  
    plot.margin          = margin(t = 15, r = 160, b = 15, l = 15),
    text            = element_text(size = 12, colour = "black"),
    axis.text       = element_text(size = 12, colour = "black")
  )

dev.off()

# Plot Bulks in the ageline to see where they lie (i.e. what does the single cell approach add?)    
# Make sure to only include autosomal, VAF>0.3 SNVs! only  

vcf_file_PRN4_bulk_lymph_node <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PRN4/batch1/snvs/batch1/PB08410_PRN4.vep_PRN4GBDLBC72/PB08410_PRN4.vep_PRN4GBDLBC72.snvs.ptato.filtered.vcf.gz"
vcf_file_P3G6_bulk_ascites <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB11197_new/intermediate/short_variants/somatic_vcfs/PB11197/231018_HMFreg2086_PB11197.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_PB11197-BLASC-BCELLBULK.SMuRF.filtered.vcf.gz"
vcf_file_P856_bulk_pleura <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB14458/intermediate/short_variants/somatic_vcfs/PB14458/231123_HMFreg2101_PB14458.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_PB14458-BLPL-BCELLBULK.SMuRF.filtered.vcf.gz"
vcf_file_P856_bulk_bone_marrow <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB14458/intermediate/short_variants/somatic_vcfs/PB14458/231123_HMFreg2101_PB14458.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_PB14458-BLBM-BCELLBULK.SMuRF.filtered.vcf.gz"
bulkWGS_vcf_files <- c(vcf_file_PRN4_bulk_lymph_node, vcf_file_P3G6_bulk_ascites, vcf_file_P856_bulk_pleura, vcf_file_P856_bulk_bone_marrow)
bulkWGS_diagnostic_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples", pattern = "*/*/*filtered.sorted.vcf.gz$", full.names = TRUE, recursive = TRUE)
bulkWGS_diagnostic_vcf_files_sub <- bulkWGS_diagnostic_vcf_files[!str_detect(bulkWGS_diagnostic_vcf_files, "PMCID211AAO.vep")] # filter out duplicate

vcf_files_bulk_all <- c(bulkWGS_vcf_files, bulkWGS_diagnostic_vcf_files_sub)
diagnostic_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_overview.csv')

pick_tumour_idx <- function(vcf, vcf_path, diag_ids) {
  cols <- colnames(vcf)
  
  if (grepl("Diagnostic_samples", vcf_path, fixed = TRUE)) {
    hit <- intersect(cols, diag_ids)
    if (length(hit) == 1) return(match(hit, cols))
    if (length(hit) > 1)
      stop("Multiple tumour columns match diagnostic IDs in ", basename(vcf_path))
    stop("No tumour column matches diagnostic IDs in ", basename(vcf_path))
  }
  
  ## non-diagnostic file: try to guess
  if (length(cols) == 1) return(1)
  
  guess <- grep("Tumor|TUMOR|_T$|-T$", cols, ignore.case = TRUE, value = TRUE)
  if (length(guess) >= 1) return(match(guess[1], cols))
  
  1  # fallback: first column
}

# ------------------------------------------------------------------ #
# helper • extract VAF vector for ONE sample column                  #
# ------------------------------------------------------------------ #
get_vaf <- function(vcf, col_idx) {
  if ("VAF" %in% names(geno(vcf))) {
    v <- geno(vcf)$VAF[, col_idx]
  } else if ("AF" %in% names(geno(vcf))) {
    v <- geno(vcf)$AF[, col_idx]
  } else if (all(c("AD", "DP") %in% names(geno(vcf)))) {
    ad <- geno(vcf)$AD[, col_idx]
    dp <- geno(vcf)$DP[, col_idx]
    v  <- ad / dp
  } else {
    stop("No VAF/AF or AD/DP fields in VCF")
  }
  as.numeric(v)
}

# autosome labels (both naming styles)
auto_plain <- as.character(1:22)
auto_chr   <- paste0("chr", 1:22)

diag_ids <- diagnostic_df$Tumor   # vector of tumour IDs to match

# ------------------------------------------------------------------ #
# main loop                                                          #
# ------------------------------------------------------------------ #
df_counts <- lapply(vcf_files_bulk_all, function(vcf_path) {
  
  vcf <- readVcf(vcf_path, genome = "hg38")
  
  col_idx <- pick_tumour_idx(vcf, vcf_path, diag_ids)
  tumour_id <- colnames(vcf)[col_idx]
  
  ## autosomes only
  keep_auto <- seqnames(vcf) %in% auto_plain | seqnames(vcf) %in% auto_chr
  vcf <- vcf[keep_auto]
  
  ## SNVs only
  is_snv <- nchar(ref(vcf)) == 1 &
    sapply(alt(vcf), \(x) length(x) == 1 && nchar(x[1]) == 1)
  vcf <- vcf[is_snv]
  
  ## VAF > 0.35
  vaf <- get_vaf(vcf, col_idx)
  vcf <- vcf[vaf > 0.35]
  
  tibble(sample = tumour_id,
         n_autosomal_SNVs_VAF_ge_0.35 = length(vcf))
}) |> bind_rows()

df_counts

