library(GenomicRanges)
library(dplyr)
library(purrr)
library(tibble)
library(readxl)
library(writexl)

# If your mts object is the flat 2×N list (GRanges, VCF, GRanges, VCF, ...):
mts <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalTimeR/Data/mts_noise_small_merged.rds")
cnvs_list <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalTimeR/Data/cnvs_list.rds")

bulk_wgs_meta <- read_xlsx("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_manuscript.xlsx")

# Re-pack into one element per sample:
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

# First: join and replace with bulk_wgs_meta IDs (as before)
cn_changes_df <- cn_changes_df %>%
  left_join(
    bulk_wgs_meta %>% select(PMC_ID, ID),
    by = c("Sample" = "PMC_ID")
  ) %>%
  mutate(Sample = if_else(!is.na(ID), ID, Sample)) %>%
  select(-ID)

# Then: apply your specific manual name fixes
cn_changes_df <- cn_changes_df %>%
  mutate(
    Sample = case_when(
      Sample == "P856_BM" ~ "P856_R2",
      Sample == "P856_PL" ~ "P856_D",
      Sample == "PMCID211AAO_2" ~ "P856_R1",
      TRUE ~ Sample
    )
  )

write_xlsx(cn_changes_df, path = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CNA/Data/cnv_data_bulk.xlsx")


