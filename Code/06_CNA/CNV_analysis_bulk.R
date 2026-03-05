################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Get CNV information for bulk WGS samples
# Author: Alexander Steemers
################################################################################

# Load libraries 
library(GenomicRanges)
library(dplyr)
library(purrr)
library(tibble)
library(readxl)
library(writexl)

# Read in input
mts <- readRDS("mts_noise_small_merged.rds") # this is generated after running MutationalTimeR.R script
cnvs_list <- readRDS("cnvs_list.rds") # this is generated after running MutationalTimeR.R script
bulk_wgs_meta <- read_xlsx("Data/Bulk_sample_manuscript.xlsx")

# Re-pack into one element per sample
mts_by_sample <- lapply(seq(1, length(mts), by = 2), function(i) {
  list(cn = mts[[i]], vcf = mts[[i + 1]])
})
names(mts_by_sample) <- names(cnvs_list)  # or another vector of sample IDs you used

cn_changes_df <- imap_dfr(mts_by_sample, function(x, samp) {
  gr <- x$cn
  if (length(gr) == 0) return(NULL)
  
  # totals (MutationTimeR stores major/minor; total = major + minor)
  tot_cn <- mcols(gr)$major_cn + mcols(gr)$minor_cn
  df <- tibble(
    Sample       = samp,
    chr          = as.character(seqnames(gr)),
    start        = start(gr),
    end          = end(gr),
    width_bp     = width(gr),
    width_Mb     = width(gr) / 1e6,
    major_cn     = mcols(gr)$major_cn,
    minor_cn     = mcols(gr)$minor_cn,
    copyNumber   = tot_cn,
    type         = mcols(gr)$type %||% NA,   # MutationTimeR event type (e.g., Mono-allelic Gain, CN-LOH)
    time         = mcols(gr)$time %||% NA,   # timing (if available)
    n_snv_mnv    = mcols(gr)$n.snv_mnv %||% NA
  ) %>%
    filter(copyNumber != 2) %>%                       # only CN changes
    filter(!chr %in% c("X","Y","chrX","chrY")) %>%    # drop sex chromosomes
    mutate(
      status = case_when(
        !is.na(type) & grepl("CN-LOH", type) ~ "cnloh",
        copyNumber > 2 ~ "gain",
        copyNumber < 2 ~ "loss",
        TRUE ~ "neutral"
      )
    )
  df
})

write_xlsx(cn_changes_df, path = "cnv_data_bulk.xlsx")


