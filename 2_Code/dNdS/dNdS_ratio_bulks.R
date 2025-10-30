################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to perform dN/dS on bulk WGS samples
# Author: Alexander Steemers
# Date: September 2025
################################################################################

# ---- Packages ----
library(VariantAnnotation)
library(dplyr)
library(purrr)
library(dndscv)
library(ggbreak)

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


# ---- Helper: normalize chromosome names to hg38 style with 'chr' ----
to_chr <- function(x) {
  x <- as.character(x)
  x <- ifelse(grepl("^chr", x, ignore.case = TRUE), x, paste0("chr", x))
  # keep only autosomes
  x[x %in% paste0("chr", 1:22)]
}

# ---- Helper: read 1 VCF -> tibble(sampleID, chr, pos, ref, alt) ----
vcf_to_mutdf <- function(vcf_path, keep_samples = NULL) {
  vcf <- VariantAnnotation::readVcf(vcf_path)
  
  # >>> this is the critical change <<<
  vcf <- VariantAnnotation::expand(vcf, row.names = FALSE)
  
  vcf_samples <- VariantAnnotation::samples(VariantAnnotation::header(vcf))
  tumor_samples <- if (is.null(keep_samples)) vcf_samples else intersect(vcf_samples, keep_samples)
  if (!length(tumor_samples)) return(NULL)
  
  rr  <- rowRanges(vcf)
  chr <- as.character(GenomeInfoDb::seqnames(rr))
  pos <- BiocGenerics::start(rr)
  ref <- as.character(VariantAnnotation::ref(vcf))
  alt <- vapply(VariantAnnotation::alt(vcf), function(a) as.character(a[1]), character(1))
  
  # autosomes only
  chr_norm <- to_chr(sub("^chr", "", chr, ignore.case = TRUE))
  keep <- !is.na(chr_norm)
  
  vcf <- vcf[keep, , drop = FALSE]
  chr <- chr_norm[keep]; pos <- pos[keep]; ref <- ref[keep]; alt <- alt[keep]
  
  purrr::map_dfr(tumor_samples, function(s) {
    gt <- if ("GT" %in% names(VariantAnnotation::geno(vcf))) VariantAnnotation::geno(vcf)$GT[, s] else rep(NA, nrow(vcf))
    k  <- !is.na(gt) & gt != "0/0" & gt != "./."
    if (!any(k)) return(NULL)
    tibble::tibble(sampleID = s, chr = chr[k], pos = as.integer(pos[k]), ref = ref[k], alt = alt[k])
  })
}

# Build mutation table from all VCFs
mut_df_list <- lapply(vcf_paths, vcf_to_mutdf)  # or vcf_to_mutdf(..., keep_samples = tumor_ids)
mut <- bind_rows(mut_df_list) %>%
  #filter(nchar(ref) == 1 & nchar(alt) == 1) %>%   # SNVs only for classic dN/dS
  distinct(sampleID, chr, pos, ref, alt, .keep_all = TRUE)

mut <- mut %>% filter(!is.na(chr))
# Quick sanity check
stopifnot(all(mut$chr %in% paste0("chr", 1:22)))
if (nrow(mut) == 0) stop("No SNVs found after filtering.")

mut <- mut %>%
  filter(chr %in% paste0("chr", 1:22))

mut$chr <- sub("^chr", "", mut$chr, ignore.case = TRUE)


# ---- Run dN/dS (hg38) ----
set.seed(1)
dndsout <- dndscv(mut, refdb = "hg38")  # uses built-in hg38 ref/annotation
sel_cv = dndsout$sel_cv
signif_genes_with_indels = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes_with_indels) = NULL
print(signif_genes_with_indels)
sel <- dndsout$sel_cv %>%
  # choose y-axis q: qglobal_cv if available, otherwise qallsubs_cv
  mutate(
    qglobal = if ("qglobal_cv" %in% names(.)) qglobal_cv else qallsubs_cv,
    # truncating effect size as the max of nonsense & splice
    wtrunc = pmax(wnon_cv, wspl_cv, na.rm = TRUE),
    # pick best (most significant) class per gene
    class_best = ifelse(qtrunc_cv <= qmis_cv, "trunc", "miss"),
    w_best  = ifelse(class_best == "trunc", wtrunc, wmis_cv),
    # plotting transforms
    log2_w  = log2(w_best),
    neglog10_q = -log10(pmax(qglobal, .Machine$double.xmin)),
    sig = qglobal < 0.10 & w_best > 1
  ) %>%
  filter(is.finite(log2_w), is.finite(neglog10_q))

# genes to label (top by q)
top_hits <- sel %>%
  arrange(qglobal) %>%
  filter(sig) %>%
  slice_head(n = 20)

ggplot(sel, aes(x = log2_w, y = neglog10_q)) +
  geom_point(aes(color = sig), alpha = 0.6) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_hline(yintercept = -log10(0.10), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = top_hits,
    aes(label = gene_name),
    size = 3, max.overlaps = 25
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "dN/dS volcano (global q, best class effect size)",
    x = "log2(dN/dS) of best class (missense vs truncating)",
    y = "-log10(global q-value)",
    color = "q < 0.10"
  )

ggplot(sel, aes(x = log2_w, y = neglog10_q)) +
  geom_point(aes(color = sig), alpha = 0.6) +
  scale_color_manual(values = c("grey70", "red")) +
  geom_hline(yintercept = -log10(0.10), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = top_hits,
    aes(label = gene_name),
    size = 3, max.overlaps = 25
  ) +
  scale_y_break(c(10, 300)) +
  labs(
    title = "dN/dS volcano (global q, best class effect size)",
    x = "log2(dN/dS) of best class (missense vs truncating)",
    y = "-log10(global q-value)",
    color = "q < 0.10"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    # remove gridlines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # bold black axis lines
    axis.line = element_line(color = "black")
  )

