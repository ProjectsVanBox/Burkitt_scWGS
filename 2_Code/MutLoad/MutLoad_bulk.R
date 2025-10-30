################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at violin plots, SNV/INDEL counts of bulk tumour samples 
# Author: Alexander Steemers
# Date: July 2025
################################################################################

# Load libraries

library(VariantAnnotation)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(readxl)
library(purrr)
library(RColorBrewer)
library(stringr)
library(patchwork)  
library(cowplot)    

# Set working directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Bulk")

# Load functions and plotting functions

source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

vcf_paths <- c(
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P3G6/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P3G6_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PRN4/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PRN4_bulk.vep.SMuRF.filtered.sorted.vcf.gz",  
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P856/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P856_D.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P856/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P856_R.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PIA9/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PIA9_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PVA9/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PVA9_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PJBU/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PJBU_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID104AAO/SMuRF/PMCID104AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID132AAL/SMuRF/PMCID132AAL.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID137AAO/SMuRF/PMCID137AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID163AAN/SMuRF/PMCID163AAN.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID321AAO/SMuRF/PMCID321AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID340AAO/SMuRF/PMCID340AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID491AAS/SMuRF/PMCID491AAS.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID509AAT/SMuRF/PMCID509AAT.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID540AAN/SMuRF/PMCID540AAN.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID610AAS/SMuRF/PMCID610AAS.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID690AAT/SMuRF/PMCID690AAT.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID821AAL/SMuRF/PMCID821AAL.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID867AAT/SMuRF/PMCID867AAT.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID967AAP/SMuRF/PMCID967AAP.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID211AAO_2/SMuRF/PMCID211AAO_2.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID458AAQ/SMuRF/PMCID458AAQ.vep.SMuRF.filtered.sorted.vcf.gz"
  )

vcf_paths <- vcf_paths[!grepl("690AAT|458AAQ", vcf_paths)] # no myc translocation found in this sample

# Function to read and filter for autosomal chromosomes
read_and_filter_vcf_all <- function(vcf_path) {
  vcf <- readVcf(vcf_path)
  chr <- as.character(seqnames(rowRanges(vcf)))
  autosomal_vcf <- vcf[!chr %in% c("X", "Y", "MT", "M", "chrX", "chrY", "chrM")]
  return(autosomal_vcf)
}

# Apply to all paths
vcf_list_all <- lapply(vcf_paths, read_and_filter_vcf_all)

# Function to read, keep autosomes, and filter to SNPs
read_and_filter_vcf_snps_only <- function(vcf_path) {
  vcf <- readVcf(vcf_path)
  
  # Keep only autosomal chromosomes
  chr <- as.character(GenomeInfoDb::seqnames(rowRanges(vcf)))
  autosomal_vcf <- vcf[!chr %in% c("X", "Y", "MT", "M", "chrX", "chrY", "chrM")]
  
  # Keep only SNPs
  snp_vcf <- autosomal_vcf[isSNV(autosomal_vcf)]
  return(snp_vcf)
}

# Apply to all paths
vcf_list_snps_only <- lapply(vcf_paths, read_and_filter_vcf_snps_only)

# Get sample names from VCF metadata
sample_names <- sapply(vcf_list_all, function(vcf) {
  samples(header(vcf))
})
sample_names <- unlist(sample_names, use.names = FALSE)

bulk_wgs_meta <- read_xlsx("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_manuscript.xlsx")

sample_names <- sample_names[sample_names %in% bulk_wgs_meta$Tumor]


# Extract VAFs from all VCFs
vaf_list <- list()

for (i in seq_along(vcf_list_all)) {
  vcf <- vcf_list_all[[i]]
  vcf_samples <- samples(header(vcf))
  
  # Identify which of these samples are tumor samples of interest
  tumor_samples <- intersect(vcf_samples, bulk_wgs_meta$Tumor)
  
  # Skip if no matching tumor samples
  if (length(tumor_samples) == 0) next
  
  # Loop through each tumor sample individually
  for (s in tumor_samples) {
    if ("AF" %in% names(geno(vcf))) {
      vaf <- geno(vcf)$AF[, s]
    } else if ("AD" %in% names(geno(vcf))) {
      ad <- geno(vcf)$AD[, s]
      vaf <- sapply(ad, function(x) {
        if (length(x) == 2 && sum(x) > 0) x[2] / sum(x) else NA
      })
    } else {
      vaf <- rep(NA, length(vcf))
    }
    
    vaf_list[[s]] <- tibble(Sample = s, VAF = vaf)
  }
}


# Combine all VAFs into one data frame
vaf_df <- bind_rows(vaf_list) %>%
  filter(!is.na(VAF), VAF != 0)

# Loop over samples and save a violin plot for each
for (s in sample_names) {
  p <- ggplot(filter(vaf_df, Sample == s), aes(x = Sample, y = VAF)) +
    coord_cartesian(ylim = c(0, 1)) +
    geom_violin(fill = "#729ECE", color = "black", scale = "width", trim = FALSE) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("VAF Distribution -", s),
      y = "Variant Allele Frequency",
      x = ""
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 5)
    )
  
  # Save to file
  ggsave(
    filename = paste0("Figures/", s, "_bulk_violin_plot.pdf"),
    plot = p,
    width = 4,
    height = 4
  )
}

# Loop over samples and save a histogram plot for each
for (s in sample_names) {
  p <- ggplot(filter(vaf_df, Sample == s), aes(x = VAF)) +
    coord_cartesian(xlim = c(0, 1)) +
    geom_histogram(
      bins = 30,              # adjust number of bins
      fill = "#729ECE",
      color = "black"
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("VAF Distribution -", s),
      y = "Count",
      x = "Variant Allele Frequency"
    ) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  
  # Save to file
  ggsave(
    filename = paste0("Figures/", s, "_bulk_histogram.pdf"),
    plot = p,
    width = 4,
    height = 4
  )
}

# Loop over samples and save a smooth histogram (density plot) for each
for (s in sample_names) {
  p <- ggplot(filter(vaf_df, Sample == s), aes(x = VAF)) +
    coord_cartesian(xlim = c(0, 1)) +
    geom_density(
      fill = "#729ECE",
      color = "black",
      alpha = 0.6
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("VAF Distribution -", s),
      y = "Density",
      x = "Variant Allele Frequency"
    ) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  
  # Save to file
  ggsave(
    filename = paste0("Figures/", s, "_bulk_density.pdf"),
    plot = p,
    width = 4,
    height = 4
  )
}

# Get Median and Mean VAF

summary_stats <- vaf_df %>%
  group_by(Sample) %>%
  summarise(
    mean_VAF   = mean(VAF, na.rm = TRUE),
    median_VAF = median(VAF, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# Now only on SNPs

sample_names_snps <- sapply(vcf_list_snps_only, function(vcf) {
  samples(header(vcf))
})
sample_names_snps <- unlist(sample_names_snps, use.names = FALSE)

sample_names_snps <- sample_names_snps[sample_names_snps %in% bulk_wgs_meta$Tumor]


# Extract VAFs from all VCFs
vaf_list_snps <- list()

for (i in seq_along(vcf_list_snps_only)) {
  vcf <- vcf_list_snps_only[[i]]
  vcf_samples <- samples(header(vcf))
  
  # Identify which of these samples are tumor samples of interest
  tumor_samples <- intersect(vcf_samples, bulk_wgs_meta$Tumor)
  
  # Skip if no matching tumor samples
  if (length(tumor_samples) == 0) next
  
  # Loop through each tumor sample individually
  for (s in tumor_samples) {
    if ("AF" %in% names(geno(vcf))) {
      vaf <- geno(vcf)$AF[, s]
    } else if ("AD" %in% names(geno(vcf))) {
      ad <- geno(vcf)$AD[, s]
      vaf <- sapply(ad, function(x) {
        if (length(x) == 2 && sum(x) > 0) x[2] / sum(x) else NA
      })
    } else {
      vaf <- rep(NA, length(vcf))
    }
    
    vaf_list_snps[[s]] <- tibble(Sample = s, VAF = vaf)
  }
}


# Combine all VAFs into one data frame
vaf_df_snps <- bind_rows(vaf_list_snps) %>%
  filter(!is.na(VAF), VAF != 0)

# Loop over samples and save a violin plot for each
for (s in sample_names) {
  p <- ggplot(filter(vaf_df_snps, Sample == s), aes(x = Sample, y = VAF)) +
    coord_cartesian(ylim = c(0, 1)) +
    geom_violin(fill = "#729ECE", color = "black", scale = "width", trim = FALSE) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("VAF Distribution -", s),
      y = "Variant Allele Frequency",
      x = ""
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 5)
    )
  
  # Save to file
  ggsave(
    filename = paste0("Figures/", s, "_bulk_violin_plot_snps.pdf"),
    plot = p,
    width = 4,
    height = 4
  )
}

# Loop over samples and save a histogram plot for each
for (s in sample_names) {
  p <- ggplot(filter(vaf_df_snps, Sample == s), aes(x = VAF)) +
    coord_cartesian(xlim = c(0, 1)) +
    geom_histogram(
      bins = 30,              # adjust number of bins
      fill = "#729ECE",
      color = "black"
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("VAF Distribution -", s),
      y = "Count",
      x = "Variant Allele Frequency"
    ) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  
  # Save to file
  ggsave(
    filename = paste0("Figures/", s, "_bulk_histogram_snps.pdf"),
    plot = p,
    width = 4,
    height = 4
  )
}

# Loop over samples and save a smooth histogram (density plot) for each
for (s in sample_names) {
  p <- ggplot(filter(vaf_df_snps, Sample == s), aes(x = VAF)) +
    coord_cartesian(xlim = c(0, 1)) +
    geom_density(
      fill = "#729ECE",
      color = "black",
      alpha = 0.6
    ) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("VAF Distribution -", s),
      y = "Density",
      x = "Variant Allele Frequency"
    ) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    )
  
  # Save to file
  ggsave(
    filename = paste0("Figures/", s, "_bulk_density_snps.pdf"),
    plot = p,
    width = 4,
    height = 4
  )
}

# Get Median and Mean VAF

summary_stats_snps <- vaf_df_snps %>%
  group_by(Sample) %>%
  summarise(
    mean_VAF   = mean(VAF, na.rm = TRUE),
    median_VAF = median(VAF, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats_snps)


# Helper: robust band from GMM; fallback to KDE if needed
band_from_gmm <- function(x, k = 2, G = 1:3, min_band = 0.05, max_band = 0.25) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(0.08)
  
  # keep inside (0,1), add tiny jitter to avoid singularities
  x <- pmin(pmax(x, 1e-4), 1 - 1e-4)
  x <- jitter(x, amount = 1e-4)
  
  fit <- tryCatch(Mclust(x, G = G, verbose = FALSE), error = function(e) NULL)
  
  if (!is.null(fit) && !is.null(fit$parameters$mean)) {
    mu <- as.numeric(fit$parameters$mean)
    k_idx <- which.max(mu)
    # SD of chosen component
    if (!is.null(fit$parameters$variance$sigma)) {
      sd_k <- sqrt(as.numeric(fit$parameters$variance$sigma[,,k_idx]))
    } else if (!is.null(fit$parameters$variance$sigmasq)) {
      sd_k <- sqrt(as.numeric(fit$parameters$variance$sigmasq[k_idx]))
    } else sd_k <- NA_real_
    
    if (is.finite(sd_k) && sd_k > 0) {
      band <- k * sd_k
      return(pmax(min_band, pmin(band, max_band)))
    }
  }
  
  # Fallback: use KDE bandwidth as spread proxy
  d <- density(x, from = 0, to = 1, na.rm = TRUE)
  bw_band <- 1.5 * d$bw
  pmax(min_band, pmin(bw_band, max_band))
}

# Get peak & band per sample with robust steps and fallback peak by KDE
gmm_peak_band <- vaf_df_snps %>%
  group_by(Sample) %>%
  group_modify(~{
    x <- .x$VAF
    x_clean <- x[is.finite(x)]
    x_clean <- pmin(pmax(x_clean, 1e-4), 1 - 1e-4)
    x_clean <- jitter(x_clean, amount = 1e-4)
    
    peak <- tryCatch({
      fit <- Mclust(x_clean, G = 1:3, verbose = FALSE)
      mu <- as.numeric(fit$parameters$mean)
      mu[which.max(mu)]
    }, error = function(e) NA_real_)
    
    if (!is.finite(peak)) {
      d <- density(x_clean, from = 0, to = 1, na.rm = TRUE)
      peak <- d$x[which.max(d$y)]
    }
    
    tibble(
      peak = peak,
      band = band_from_gmm(x_clean, k = 2)
    )
  }) %>%
  ungroup()

# Recompute flags & clonal set
clonal_flagged <- vaf_df_snps %>%
  inner_join(gmm_peak_band, by = "Sample") %>%
  mutate(clonal = is.finite(peak) & abs(VAF - peak) <= band)

clonal_vafs <- clonal_flagged %>% filter(clonal) %>% dplyr::select(-clonal)

# --- Build per-sample plotting data with peak and band ------------------------
vaf_plot_df <- vaf_df_snps %>%
  dplyr::inner_join(gmm_peak_band, by = "Sample") %>%
  dplyr::mutate(
    band_lo = pmax(0, peak - band),
    band_hi = pmin(1, peak + band)
  )

p_vaf_facets <- ggplot(vaf_plot_df, aes(x = VAF)) +
  # clonal band as background ribbon (per-sample)
  geom_rect(
    data = vaf_plot_df %>% dplyr::distinct(Sample, band_lo, band_hi),
    aes(xmin = band_lo, xmax = band_hi, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, alpha = 0.15, fill = "firebrick"
  ) +
  geom_histogram(bins = 60, fill = "grey80", color = "grey60") +
  geom_density(linewidth = 0.6) +
  # peak line (per-sample)
  geom_vline(
    data = vaf_plot_df %>% dplyr::distinct(Sample, peak),
    aes(xintercept = peak),
    color = "firebrick", linewidth = 0.8
  ) +
  facet_wrap(~ Sample, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "VAF distributions per sample with clonal band",
    x = "Variant Allele Fraction (VAF)",
    y = "Count / Density"
  ) +
  theme_bw(base_size = 11)

p_vaf_facets

clonal_counts <- clonal_flagged %>%
  mutate(clonal_i = as.integer(clonal)) %>%
  dplyr::count(Sample, wt = clonal_i, name = "n_clonal")

write.csv(clonal_counts, "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Bulk/Data/clonal_counts_per_sample.csv", row.names = FALSE)

# Exclusions (as before)
excluded <- c("PMCID90AAT", "PMCID458AAQ")

# Clean inputs
bulk_wgs_meta_clean <- bulk_wgs_meta %>%
  mutate(across(everything(), ~ if (is.character(.x) || is.factor(.x)) str_trim(as.character(.x)) else .x)) %>%
  filter(!PMC_ID %in% excluded)

gmm_peak_band <- gmm_peak_band %>%
  mutate(Sample = str_trim(as.character(Sample))) %>%
  filter(!Sample %in% excluded)

vaf_df_snps <- vaf_df_snps %>%
  mutate(Sample = str_trim(as.character(Sample))) %>%
  filter(!Sample %in% excluded)

all_samples <- gmm_peak_band %>% distinct(Sample) %>% pull(Sample)

# --- 1) Primary lookup: PMC_ID -> ID (collapsing duplicates) ---
id_lookup_pmc <- bulk_wgs_meta_clean %>%
  filter(!is.na(PMC_ID), PMC_ID != "", !is.na(ID), ID != "") %>%
  group_by(PMC_ID) %>%
  summarise(ID = paste(unique(ID), collapse = " / "), .groups = "drop") %>%
  deframe()

# --- 2) Heuristic: find the meta column that best matches your Sample values ---
# Consider only character-like columns
candidate_cols <- names(bulk_wgs_meta_clean)[
  vapply(bulk_wgs_meta_clean, function(x) is.character(x) || is.factor(x), logical(1))
]

coverage <- map_dfr(candidate_cols, function(col) {
  vals <- unique(str_trim(as.character(bulk_wgs_meta_clean[[col]])))
  tibble(column = col, coverage = mean(all_samples %in% vals, na.rm = TRUE))
}) %>%
  arrange(desc(coverage))

best_col <- coverage %>%
  filter(coverage > 0) %>%
  slice(1) %>%
  pull(column)

if (length(best_col) == 0) {
  message("No matching column found in bulk_wgs_meta for Sample values; using PMC_ID-only mapping.")
  # Build labels using PMC_ID-only
  sample_labels <- tibble(Sample = all_samples) %>%
    mutate(display = unname(id_lookup_pmc[Sample]),
           display = ifelse(is.na(display) | display == "", Sample, display))
} else {
  message("Using `", best_col, "` to map Sample -> ID (coverage: ",
          round(coverage$coverage[coverage$column == best_col] * 100, 1), "%)")
  
  # Build secondary lookup: best_col -> ID
  id_lookup_best <- bulk_wgs_meta_clean %>%
    filter(!is.na(.data[[best_col]]), .data[[best_col]] != "", !is.na(ID), ID != "") %>%
    group_by(.data[[best_col]]) %>%
    summarise(ID = paste(unique(ID), collapse = " / "), .groups = "drop") %>%
    { setNames(.$ID, .[[best_col]]) }
  
  # Labels: try PMC_ID -> ID first, then best_col -> ID, else fallback to Sample
  sample_labels <- tibble(Sample = all_samples) %>%
    mutate(
      lab_pmc  = unname(id_lookup_pmc[Sample]),
      lab_alt  = unname(id_lookup_best[Sample]),
      display  = dplyr::coalesce(ifelse(lab_pmc == "", NA, lab_pmc),
                                 ifelse(lab_alt == "", NA, lab_alt),
                                 Sample)
    ) %>%
    select(Sample, display)
}

# ---- Use sample_labels$display for titles, as before ----
plot_vaf_density_with_band_label <- function(sample_id, display_label, vaf_df_snps, gmm_peak_band) {
  df <- dplyr::filter(vaf_df_snps, Sample == sample_id)
  pb <- dplyr::filter(gmm_peak_band, Sample == sample_id)
  if (nrow(pb) == 0 || !is.finite(pb$peak[1]) || !is.finite(pb$band[1])) return(NULL)
  
  peak <- pb$peak[1]; band <- pb$band[1]
  xmin_band <- max(0, peak - band); xmax_band <- min(1, peak + band)
  
  ggplot(df, aes(x = VAF)) +
    annotate("rect", xmin = xmin_band, xmax = xmax_band, ymin = -Inf, ymax = Inf,
             alpha = 0.2, fill = "#E69F00") +
    geom_density(fill = "#729ECE", color = "black", alpha = 0.6) +
    geom_vline(xintercept = peak, linetype = "dashed", size = 0.7) +
    coord_cartesian(xlim = c(0, 1)) +
    theme_minimal(base_size = 10) +
    labs(
      title = paste0(display_label),
      x = "Variant Allele Frequency", y = "Density"
    )
}

plots <- lapply(seq_len(nrow(sample_labels)), function(i) {
  sid <- sample_labels$Sample[i]
  lab <- sample_labels$display[i]
  plot_vaf_density_with_band_label(sid, lab, vaf_df_snps, gmm_peak_band)
})
plots <- Filter(Negate(is.null), plots)

# (If you still see old IDs, print diagnostics)
print(coverage)
print(head(sample_labels, 10))

# 1) Build plots using the mapped display labels
plots <- map2(
  sample_labels$Sample,
  sample_labels$display,
  ~ plot_vaf_density_with_band_label(.x, .y, vaf_df_snps, gmm_peak_band)
)

# Drop NULLs (samples without peak/band)
plots <- keep(plots, ~ !is.null(.x))

# 2) Ensure exactly 3 cols × 7 rows (21 panels) by padding with blank panels if needed
need <- 21 - length(plots)
if (need > 0) {
  blank <- ggplot() + theme_void() + theme(plot.margin = margin(2,2,2,2))
  plots <- c(plots, rep(list(blank), need))
}

# 3) Arrange into grid and add legend space
grid <- wrap_plots(plots, ncol = 3)

final_plot <- grid
# 4) Save to A4 (portrait)
ggsave(
  filename = "Figures/Clonal_selection_bulk.pdf",
  plot = final_plot,
  width = 8.27, height = 8, units = "in"
)

#Get variant types (SNP or INDEL) per sample
mutation_summary <- map_dfr(seq_along(vcf_list_all), function(i) {
  vcf <- vcf_list_all[[i]]
  vcf_samples <- samples(header(vcf))
  tumor_samples <- intersect(vcf_samples, bulk_wgs_meta$Tumor)
  if (length(tumor_samples) == 0) return(NULL)
  
  ref <- as.character(ref(vcf))
  alt <- sapply(alt(vcf), function(x) as.character(x[[1]]))
  variant_type <- ifelse(nchar(ref) == 1 & nchar(alt) == 1, "SNP", "INDEL")
  
  # Loop through tumor samples and only keep variants with genotype calls
  map_dfr(tumor_samples, function(s) {
    # Genotype field for this sample
    gt <- geno(vcf)$GT[, s]
    # Keep only variants with a call (e.g. not "./.")
    keep <- !is.na(gt) & gt != "0/0" & gt != "./."
    
    tibble(
      Sample = s,
      Variant_Type = variant_type[keep]
    )
  })
})

# Now proceed with counting
mutation_counts <- mutation_summary %>%
  dplyr::count(Sample, Variant_Type) %>%
  tidyr::pivot_wider(names_from = Variant_Type, values_from = n, values_fill = 0)

mutation_stats <- mutation_counts %>%
  summarise(
    SNP_median = median(SNP, na.rm = TRUE),
    SNP_range = paste0(min(SNP, na.rm = TRUE), "–", max(SNP, na.rm = TRUE)),
    INDEL_median = median(INDEL, na.rm = TRUE),
    INDEL_range = paste0(min(INDEL, na.rm = TRUE), "–", max(INDEL, na.rm = TRUE))
  )

print(mutation_stats)

sum(mutation_counts$SNP) + sum(mutation_counts$INDEL)

mutation_counts <- mutation_counts %>%
  left_join(bulk_wgs_meta, by = c("Sample" = "Tumor"))

mutation_counts <- mutation_counts %>%
  filter(Timepoint != "Relapse") # take only diagnosis

mutation_counts <- mutation_counts %>%
  filter(Sample != "PMLBM000KOD") # because it's an outlier!

ggplot(mutation_counts, aes(x = Age_at_diagnosis, y = SNP)) +
  geom_point(color = "#1f78b4", size = 2) +
  theme_CHemALL() + 
  scale_x_continuous(limits = c(0, 18), breaks = seq(0, 18, by = 2))+
  geom_smooth(method = "lm", color = "black", linetype = "solid", se = FALSE, fullrange = TRUE) +
  labs(
    title = "SNV Count vs Age at Diagnosis",
    x = "Age at Diagnosis",
    y = "SNV Count"
  ) 

ggplot(mutation_counts, aes(x = Age_at_diagnosis, y = INDEL)) +
  geom_point(color = "#1f78b4", size = 2) +
  theme_CHemALL() + 
  scale_x_continuous(limits = c(0, 18), breaks = seq(0, 18, by = 2))+
  geom_smooth(method = "lm", color = "black", linetype = "solid", se = FALSE, fullrange = TRUE) +
  labs(
    title = "INDEL Count vs Age at Diagnosis",
    x = "Age at Diagnosis",
    y = "INDEL Count"
  ) 
# Linear model for SNVs
lm_snp <- lm(SNP ~ Age_at_diagnosis, data = mutation_counts)
summary(lm_snp)

# Linear model for INDELs
lm_indel <- lm(INDEL ~ Age_at_diagnosis, data = mutation_counts)
summary(lm_indel)

# Extract p-value and R2 for SNPs
snp_summary <- summary(lm_snp)
snp_pval <- signif(snp_summary$coefficients[2, 4], 3)
snp_r2 <- signif(snp_summary$r.squared, 3)

# Extract p-value and R2 for INDELs
indel_summary <- summary(lm_indel)
indel_pval <- signif(indel_summary$coefficients[2, 4], 3)
indel_r2 <- signif(indel_summary$r.squared, 3)

pdf("Figures/SNV_age_corr_diagnosis_samples_only.pdf") 
ggplot(mutation_counts, aes(x = Age_at_diagnosis, y = SNP)) +
  geom_point(color = "#1f78b4", size = 2) +
  theme_CHemALL() + 
  scale_x_continuous(limits = c(0, 18), breaks = seq(0, 18, by = 2)) +
  geom_smooth(method = "lm", color = "black", linetype = "solid", se = FALSE, fullrange = TRUE) +
  labs(
    x = "Age at Diagnosis",
    y = "Somatic SNV Count"
  ) +
  annotate("text", x = 1, y = max(mutation_counts$SNP, na.rm = TRUE), 
           label = paste0("p = ", snp_pval, "\nR² = ", snp_r2), 
           hjust = 0, vjust = 1, size = 5) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
dev.off()

pdf("Figures/INDEL_age_corr_diagnosis_samples_only.pdf") 
ggplot(mutation_counts, aes(x = Age_at_diagnosis, y = INDEL)) +
  geom_point(color = "#1f78b4", size = 2) +
  theme_CHemALL() + 
  scale_x_continuous(limits = c(0, 18), breaks = seq(0, 18, by = 2)) +
  geom_smooth(method = "lm", color = "black", linetype = "solid", se = FALSE, fullrange = TRUE) +
  labs(
    x = "Age at Diagnosis",
    y = "Somatic INDEL Count"
  ) +
  annotate("text", x = 1, y = max(mutation_counts$INDEL, na.rm = TRUE), 
           label = paste0("p = ", indel_pval, "\nR² = ", indel_r2), 
           hjust = 0, vjust = 1, size = 5) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
dev.off()

# Run t-test for SNP - sex corr
snp_sex_ttest <- t.test(SNP ~ Sex, data = mutation_counts)
snp_pval <- signif(snp_sex_ttest$p.value, 3)  # round nicely

# Plot with p-value annotation
pdf("Figures/SNP_sex_corr_diagnosis_samples_only.pdf") 
ggplot(mutation_counts, aes(x = Sex, y = SNP, fill = Sex)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_CHemALL() + 
  labs(x = "Sex", y = "Somatic SNP Count") +
  annotate("text",
           x = 1.5,  # midpoint between the two boxes
           y = max(mutation_counts$SNP, na.rm = TRUE) * 1.05, # a bit above top
           label = paste0("t-test p = ", snp_pval),
           size = 5,
           hjust = 0.5) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
)
dev.off()

# Run t-test for INDEL - sex corr
indel_sex_ttest <- t.test(INDEL ~ Sex, data = mutation_counts)
indel_pval <- signif(indel_sex_ttest$p.value, 3)  # round nicely

# Plot with p-value annotation
pdf("Figures/INDEL_sex_corr_diagnosis_samples_only.pdf") 
ggplot(mutation_counts, aes(x = Sex, y = INDEL, fill = Sex)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  theme_CHemALL() + 
  labs(x = "Sex", y = "Somatic INDEL Count") +
  annotate("text",
           x = 1.5,  # midpoint between the two boxes
           y = max(mutation_counts$INDEL, na.rm = TRUE) * 1.05, # a bit above top
           label = paste0("t-test p = ", indel_pval),
           size = 5,
           hjust = 0.5) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
dev.off()


#  For both SNVs and INDELS - Slight positive trend, but not statistically significant (p-value) and the variation (R-squared) is large
# i.e. Age does not strongly predict mutation load (INDELs or SNVs) in this dataset.
# Also sex doesn't affect mutational load
# Sample size may be limiting power to detect an effect

# Empty list to collect data

# Steps: Take only PASS variants --> only SNVs --> only autosomal --> VAF >= 0.15

all_patients_data <- list()

get_ref <- function(x) if (!is.null(x) && length(x) >= 1) as.numeric(x[1]) else NA
get_alt <- function(x) if (!is.null(x) && length(x) >= 2) as.numeric(x[2]) else NA

for (Patient in Patients) {
  
  #vcf_file_unzipped <- paste0(vcf_base, Patient, vcf_suffix, Patient, ".vep.vcf")
  #vcf_file_gzipped <- paste0(vcf_base, Patient, vcf_suffix, Patient, ".vep.vcf.gz")
  
  vcf_file_unzipped <- paste0(vcf_base, Patient, vcf_suffix, Patient, "_bulk.vep.SMuRF.filtered.sorted.VAF03.vcf")
  vcf_file_gzipped <- paste0(vcf_base, Patient, vcf_suffix, Patient, "_bulk.vep.SMuRF.filtered.sorted.VAF03.vcf.gz")
  
  if (file.exists(vcf_file_unzipped)) {
    vcf_file <- vcf_file_unzipped
  } else if (file.exists(vcf_file_gzipped)) {
    vcf_file <- vcf_file_gzipped
  } else {
    cat("No VCF file found for", Patient, "\n")
    next
  }
  
  vcf <- readVcf(vcf_file, "hg38")
  
  # Keep only PASS variants
  vcf <- vcf[as.vector(unlist(fixed(vcf)$FILTER)) == "PASS"]
  
  # Keep only SNVs (remove indels)
  ref_allele <- as.character(ref(vcf))
  alt_list <- alt(vcf)
  
  is_snv <- mapply(function(ref, alt) {
    nchar(ref) == 1 && all(nchar(as.character(alt)) == 1)
  }, ref_allele, alt_list)
  
  #vcf <- vcf[is_snv]
  
  num_pass <- length(vcf)
  cat("Patient:", Patient, "- PASS variants:", num_pass, "\n")
  
  # Autosomes only
  vcf <- vcf[seqnames(rowRanges(vcf)) %in% as.character(1:22)]
  
  # REF / ALT depths
  AD <- geno(vcf)$AD
  ref_depth <- matrix(sapply(AD, get_ref), nrow = dim(AD)[1], ncol = dim(AD)[2])
  alt_depth <- matrix(sapply(AD, get_alt), nrow = dim(AD)[1], ncol = dim(AD)[2])
  
  vaf <- alt_depth / (ref_depth + alt_depth)
  vaf[vaf < 0.15 | is.na(vaf)] <- NA
  
  # Variant names
  variant_ids <- paste0(seqnames(rowRanges(vcf)), ":", start(rowRanges(vcf)))
  
  # Prepare dataframe
  vaf_df <- as.data.frame(vaf)
  colnames(vaf_df) <- colnames(vcf)
  vaf_df$Variant <- variant_ids
  
  # Long format
  vaf_long <- pivot_longer(vaf_df, cols = -Variant, names_to = "Sample", values_to = "VAF") %>%
    filter(!is.na(VAF)) %>%
    mutate(Donor = Patient)
  
  all_patients_data[[Patient]] <- vaf_long
  
  cat("Done with Patient:", Patient, "\n")
}

# Combine all patients

combined_df <- do.call(rbind, all_patients_data)
PRN4_P1B11_df <- combined_df[combined_df$Sample == "PB08410-BLLN-BCELLP1B11", ]
P3G6_P1B4_df <- combined_df[combined_df$Sample == "PB11197-BLASC-BCELLP1B4", ]
P3G6_msc_df <- combined_df[combined_df$Sample == "PB11197-BLBM-MSCBULK", ]

# Define bulk samples

bulk_WGS <- c("PB11197-BLASC-BCELLBULK", "PRN4GBDLBC72", "PB14458-BLPL-BCELLBULK","PB14458-BLBM-BCELLBULK", "PIA9GBDABC78", "PVA9GBDABC78","PJBUGBDABC82")

# Filter out MSCs

combined_df <- combined_df %>%
  filter(!grepl("MS", Sample))

# Separate bulk and single-cell

bulk_df <- combined_df %>% filter(Sample %in% bulk_WGS)
single_cell_df <- combined_df %>% filter(!Sample %in% bulk_WGS)

# Violin plots with cutoffs for bulk 30X samples

for (donor in Patients) {
  
  donor_df <- bulk_df %>% filter(Donor == donor)
  
  # Decide intercept
  if ("PRN4GBDLBC72" %in% donor_df$Sample) {
    intercept <- 0.25
  } else {
    intercept <- 0.3
  }
  
  p <- ggplot(donor_df, aes(x = Sample, y = VAF, fill = Sample)) +
    geom_violin(trim = FALSE, scale = "width", fill = "steelblue") +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white", color = "black") +
    theme_bw() +
    geom_hline(yintercept = intercept, colour = "red") +
    scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
  theme(axis.text.x = element_text(),
          strip.text = element_text(size = 10),
          strip.background = element_rect(fill = "grey90", color = NA)) +
    labs(title = paste0("Patient ", donor, " - Bulk sample(s)"),
         y = "VAF", x = "Sample") +
    guides(fill = "none")
  
  print(p)
  
  ggsave(paste0("Figures/VAF_distribution_bulk_", donor, ".pdf"), plot = p, width = 7, height = 4)
}

# Apply bulk-specific VAF cutoffs (determined by violin plots)

bulk_df <- bulk_df %>%
  rowwise() %>%
  mutate(VAF = ifelse(
    Sample == "PRN4GBDLBC72" & VAF < 0.25, NA,
    ifelse(Sample != "PRN4GBDLBC72" & VAF < 0.3, NA, VAF)
  )) %>%
  filter(!is.na(VAF)) %>%
  ungroup()

# Calculate purity of bulks (determine which are uncontaminated and which are contaminated)

bulk_purity <- bulk_df %>%
  group_by(Sample) %>%
  summarise(median_VAF = median(VAF, na.rm = TRUE)) %>%
  arrange(desc(median_VAF))

bulk_purity$Purity <- bulk_purity$median_VAF*2 # PRN4 LN contaminated --> adjust VAF cutoff
print(bulk_purity)

# Violin plots with cutoffs for single cells 15X samples

single_cell_df <- single_cell_df %>%
  mutate(Sample = factor(Sample),
         Donor = factor(Donor))

unique_donors <- unique(single_cell_df$Donor)

for (donor in unique_donors) {
  
  donor_df <- single_cell_df %>% 
    filter(Donor == donor)
  
  # Order samples by median VAF (high to low)
  sample_order <- donor_df %>%
    group_by(Sample) %>%
    summarise(median_VAF = median(VAF, na.rm = TRUE)) %>%
    arrange(desc(median_VAF)) %>%
    pull(Sample)
  
  donor_df <- donor_df %>%
    mutate(Sample = factor(Sample, levels = sample_order))
  
  p <- ggplot(donor_df, aes(x = Sample, y = VAF)) +
    geom_violin(fill = "steelblue", color = "black", trim = FALSE, scale = "width") +
    stat_summary(fun = median, geom = "point", shape = 23, size = 2, fill = "white", color = "black") +
    theme_bw() +
    geom_hline(yintercept = 0.15, color = "red") +
    scale_y_continuous(limits = c(0, 1), oob = scales::squish) +
    theme(
      axis.text.x = element_text(size = 4,angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    ) +
    labs(
      title = paste0("Patient ", donor, " - Single-cell samples"),
      y = "VAF",
      x = "Sample"
    )
  
  print(p)
  
  ggsave(paste0("Figures/VAF_distribution_single_cell_", donor, ".pdf"), plot = p, width = 7, height = 4)
  
}

# Save RDS files

saveRDS(bulk_df, "../MutLoad/Data/PASS_autosomal_VAFcutoff_bulk_samples.RDS")
saveRDS(single_cell_df, "../MutLoad/Data/PASS_autosomal_VAFcutoff_single_cell_samples.RDS")

# Optional read RDS here

bulk_df <- readRDS(file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Data/PASS_autosomal_VAFcutoff_bulk_samples.RDS")
single_cell_df <- readRDS(file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Data/PASS_autosomal_VAFcutoff_single_cell_samples.RDS")
          
# Remove below curve samples and determine VAF cutoff using MAD outlier detection

below_curve_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv")
low_call_frac_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_frac_samples.csv")

single_cell_df_filtered <- single_cell_df %>%
  filter(!Sample %in% below_curve_df$Sample_name,
         !Sample %in% low_call_frac_df)

# Parameters

num_bins <- 10
epsilon <- 1e-6  # to avoid zero-probability issues

# Decide here if using filtered data (coverage quality) or unfiltered data (all samples)

# 1. Bin VAFs into histogram for each sample
binned_df <- single_cell_df_filtered %>%
  mutate(bin = cut(VAF, breaks = seq(0, 1, length.out = num_bins + 1), include.lowest = TRUE)) %>%
  group_by(Sample, bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(prob = (count + epsilon) / sum(count + epsilon)) %>%  # normalize
  ungroup()

# 2. Pivot to wide format: one row per sample
wide_probs <- binned_df %>%
  dplyr::select(Sample, bin, prob) %>%
  pivot_wider(names_from = bin, values_from = prob, values_fill = list(prob = epsilon)) %>%
  column_to_rownames("Sample")

# 3. Compute reference distribution (e.g., median across samples)
ref_dist <- apply(wide_probs, 2, median)

# 4. Compute TVD for each sample
tvd <- function(p, q) {
  0.5 * sum(abs(p - q))
}
tvd_values <- apply(wide_probs, 1, function(p) tvd(p, ref_dist))

# 5. Output: samples ranked by TVD
tvd_df <- data.frame(Sample = names(tvd_values), TVD = tvd_values) %>%
  arrange(desc(TVD))

# 6. Plot TVD scores
mad_val <- mad(tvd_df$TVD)
median_val <- median(tvd_df$TVD)

tvd_df <- tvd_df %>%
  mutate(Flagged = TVD > (median_val + 2.5 * mad_val)) # https://www.sciencedirect.com/science/article/pii/S0022103113000668?via%3Dihub

ggplot(tvd_df, aes(x = TVD, y = reorder(Sample, TVD), 
                   fill = Flagged)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Total Variation Distance of VAF Distributions",
       x = "TVD", y = "Sample_name") +
  scale_fill_manual(values = c('#54BFB7','grey')) + theme_CHemALL() +
  theme(text =  element_text(size =  7, color = 'black'),
        axis.text = element_text(angle = 90, size = 5, colour = "black"),
        axis.text.y = element_text(angle = 90, hjust = 0.5))

# Annotate the other plot with this
plot_df3b <- merge(single_cell_df, tvd_df)

# Remove flagged with higher than median VAF
median(plot_df3b$VAF)
median_df <- plot_df3b %>% group_by(Sample) %>% summarise(med = median(VAF))

plot_df3b[plot_df3b$Sample %in% median_df[median_df$med > median(plot_df3b$VAF),]$Sample, 'Flagged'] <- FALSE

# Rename that column
plot_df3b$VAFfilter <- 'Pass'
plot_df3b[plot_df3b$Flagged,]$VAFfilter <- 'Fail'

unique_donors <- unique(plot_df3b$Donor)

for (donor in unique_donors) {
  
  donor_df <- plot_df3b %>% 
    filter(Donor == donor)
  
  # Calculate median VAF per sample
  sample_medians <- donor_df %>% 
    group_by(Sample) %>% 
    summarise(median_vaf = median(VAF, na.rm = TRUE)) %>% 
    arrange(desc(median_vaf))
  
  # Reorder Sample factor based on median VAF
  donor_df$Sample <- factor(donor_df$Sample, levels = sample_medians$Sample)
  
  # Plot
  p <- ggplot(data = donor_df,
              aes(x = Sample,
                  y = VAF,
                  fill = VAFfilter)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.8) +
    ggtitle(paste0('VAF distribution per single-cell sample - ', donor)) +
    scale_fill_manual(values = c('Fail' = 'grey', 'Pass' = '#54BFB7')) +
    theme_CHemALL() +
    ggTextAxisRotate() +
    theme(text = element_text(size = 7, color = 'black'),
          axis.text = element_text(size = 5, colour = "black"))
  
  print(p)
}

fail_samples <- plot_df3b %>%
  filter(VAFfilter == "Fail") %>%
  pull(Sample) %>%
  unique() %>%
  as.character()

blacklist_samples <- unique(c(filtered_out, fail_samples, below_curve_df$Sample_name))

P856_blacklist <- blacklist_samples[grepl("^P856|^PB14458", blacklist_samples)]
P3G6_blacklist <- blacklist_samples[grepl("^P3G6|^PB11197", blacklist_samples)]
PRN4_blacklist <- blacklist_samples[grepl("^PRN4|^PB08410", blacklist_samples)]
PIA9_blacklist <- blacklist_samples[grepl("^PIA9", blacklist_samples)]
PVA9_blacklist <- blacklist_samples[grepl("^PVA9", blacklist_samples)]
PJBU_blacklist <- blacklist_samples[grepl("^PJBU", blacklist_samples)]

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx')
colnames(input_df)
input_df_sc <- input_df[input_df$ResolveDNA_version %in% c("v1", "v2", "v2.0"), ]

perc_removed <- (length(blacklist_samples))/length(input_df_sc$Sample_name)*100
print(perc_removed)

abs_num_kept <- length(input_df_sc$Sample_name) - length(blacklist_samples)

input_df_sub <- input_df_sc[!(input_df_sc$Sample_name %in% blacklist_samples), ]

input_df_sub %>%
  group_by(Novogene_ID) %>%
  summarise(Number_of_sc = n())
