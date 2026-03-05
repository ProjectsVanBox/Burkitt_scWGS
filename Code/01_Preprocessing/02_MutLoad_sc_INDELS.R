################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Script to get filtered SNVs (total + autosomal) from single cell WGS samples 
# Author: Alexander Steemers
################################################################################

# Load libraries

library(reshape2)
library(ggplot2)
library(tidyverse)
library(VariantAnnotation)
library(readxl)
library(BSgenome)
library(GenomicRanges)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)

# Load functions and plotting functions

source('Functions/GenericFunctions.R')
source('Functions/theme_burkitt.R')

# Load metadata

input_df <-  read_excel('Data/Sample_overview.xlsx')
input_df_sc <- input_df[input_df$ResolveDNA_version %in% c("v1", "v2", "v2.0"), ]
input_df_sub <- input_df_sc[!is.na(input_df_sc$Callable_fraction) & !is.na(input_df_sc$Mean_coverage),]

# Get ptato filtered vcf files (these should be downloaded from Mendeley DOI: 10.17632/phk5jhhm7d.1)
# Find out which VFCs to download by checking the file names below
# They should be 332 single cells in total

single_cell_file_names <- readLines("Data/single_cell_ptato_filtered_file_names_INDELS.txt")

#all_filtered_vcfs <- "~/path/to/vcf/files.vcf"

# MinimalVAF cut-off

MinimalVAF <- 0.15

# Loop over all samples and make a list 

INDELS_raw <- list()

for (Sample in input_df_sub$Sample_name) {
  message("→ processing ", Sample)
  
  filtered_vcf_path   <- all_filtered_vcfs[grepl(Sample, all_filtered_vcfs, ignore.case = TRUE)][1]
  
  print(filtered_vcf_path) # to check if the right VCF file was used
  
  vcf <- readVcf(filtered_vcf_path)
  
  # Read VAF
  if (!"VAF" %in% names(geno(vcf)))
    stop("`geno(vcf)$VAF` not present in ", basename(filtered_vcf_path))
  
  vaf <- geno(vcf)$VAF
  if (length(dim(vaf)) == 2L)       
    vaf <- vaf[, 1, drop = TRUE]
  
  # Add FILTER logic
  rr <- rowRanges(vcf)
  mcols(rr)$Chromosome <- as.character(seqnames(rr))
  mcols(rr)$VAF    <- vaf
  mcols(rr)$FILTER <- ifelse(vaf >= MinimalVAF, "PASS", "FAIL_VAF")
  
  INDELS_raw[[Sample]]  <- rr
}

# Save both autosomal and sex chromosome SNVs

#saveRDS(INDELS_raw, file = "total_indel_variants_ResolveDNA_VAF015.RDS") # make sure to change minimal VAF

# Read RDS files 

INDELS_raw_015 <- readRDS(file = "Data/total_indel_variants_ResolveDNA_VAF015.RDS")

# Only autosomal
auto_chrs <- as.character(1:22)

INDELS_raw_015_autosomal <-  lapply(INDELS_raw_015, function(gr) {
  keep <- seqnames(gr) %in% auto_chrs
  gr2  <- gr[keep]                        
  keepSeqlevels(gr2, auto_chrs, pruning.mode = "coarse")
}
)

# Remove all variants that are UNCALLABLE and/or with a PTAprobs < PTAprobsCutoff
INDELS_PASS_015_autosomal <- lapply(INDELS_raw_015_autosomal, function(x) x[which(x$FILTER =="PASS"),])
INDELS_FAIL_VAF_015_autosomal <- lapply(INDELS_raw_015_autosomal, function(x) x[which(x$FILTER =="FAIL_VAF"),])

## save R objects 
#saveRDS(INDELS_PASS_015_autosomal, file = "autosomal_INDEL_PASS_variants_VAF015.RDS")
#saveRDS(INDELS_FAIL_VAF_015_autosomal, file = "autosomal_INDEL_FAILVAF_variants_VAF015.RDS")

# Get lengths for each SBS label per sample

get_lengths <- function(lst) vapply(lst, length, integer(1))

counts <- tibble(
  Sample    = names(INDELS_PASS_015_autosomal),                     
  PASS      = get_lengths(INDELS_PASS_015_autosomal),
  FAIL      = get_lengths(INDELS_FAIL_VAF_015_autosomal)
) |>
  pivot_longer(-Sample, names_to = "Category", values_to = "Count")

counts <- counts %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# Plot stacked barplot to see how much is filtered out with Minimal cutoff of 0.2

ggplot(counts, aes(x = Sample, y = Count, fill = Category)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = c(PASS = "forestgreen", FAIL = "firebrick")) +
  labs(x = "Sample", y = "Number of autosomal SNVs", fill = "Filter status") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1))

#ggsave(paste0("Mut_load_INDEL_filtering_stacked_barplot_VAF015.pdf"))
