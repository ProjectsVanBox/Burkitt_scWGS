################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Script to perform dN/dS on bulk WGS samples
# Author: Alexander Steemers
################################################################################

#Load libraries
library(VariantAnnotation)
library(dplyr)
library(purrr)
library(dndscv)
library(ggbreak)
library(ggplot2)

# For vcf_files take all "*.vep.SMuRF.filtered.sorted.vcf.gz" files found from Mendeley doi: 10.17632/phk5jhhm7d.1 in the Bulk samples and Diagnostic samples folder

# Function to normalize chromosome names to hg38 style with 'chr'
to_chr <- function(x) {
  x <- as.character(x)
  x <- ifelse(grepl("^chr", x, ignore.case = TRUE), x, paste0("chr", x))
  # keep only autosomes
  x[x %in% paste0("chr", 1:22)]
}

# Function to read VCF as a mutation df
vcf_to_mutdf <- function(vcf_path, keep_samples = NULL) {
  vcf <- VariantAnnotation::readVcf(vcf_path)
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
    k  <- !is.na(gt) & gt != "0/0" & gt != "0|0" & gt != "./."
    if (!any(k)) return(NULL)
    tibble::tibble(sampleID = s, chr = chr[k], pos = as.integer(pos[k]), ref = ref[k], alt = alt[k])
  })
}

# Build mutation table from all VCFs
mut_df_list <- lapply(vcf_paths, vcf_to_mutdf)
mut <- bind_rows(mut_df_list) %>%
  #filter(nchar(ref) == 1 & nchar(alt) == 1) %>%   # SNVs only for classic dN/dS
  distinct(sampleID, chr, pos, ref, alt, .keep_all = TRUE)

mut <- mut %>% filter(!is.na(chr))

mut$chr <- sub("^chr", "", mut$chr, ignore.case = TRUE)

# Run dN/dS
set.seed(1)
dndsout <- dndscv(mut, refdb = "hg38") 
sel_cv = dndsout$sel_cv
signif_genes_with_indels = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes_with_indels) = NULL
print(signif_genes_with_indels)
sel <- dndsout$sel_cv %>%
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

# Plot genes under positive selection
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

# Plot again, but we the y-axis now placing the outlier (TP53) closer to the other genes
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