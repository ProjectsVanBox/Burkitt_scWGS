################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Cohort CNV frequency circos 
# Author: Alexander Steemers
################################################################################

# Load libraries 
library(data.table)
library(GenomicRanges)
library(circlize)

# Define functions
process_one <- function(x, merge_gap, si) {
  # Accepts NULL, data.frame/data.table, or GRanges
  if (is.null(x)) return(make_empty(si))
  if (is.data.frame(x)) {
    if (nrow(x) == 0L) return(make_empty(si))
    x <- to_gr_from_dt(x)
  }
  if (length(x) == 0L) return(make_empty(si))
  x <- reduce(x, min.gapwidth = merge_gap)
  keep <- intersect(seqlevels(x), seqlevels(si))
  if (length(keep)) x <- keepSeqlevels(x, keep, pruning.mode = "coarse") else x <- GRanges()
  seqlevels(x) <- seqlevels(si)
  suppressWarnings(seqinfo(x) <- si)
  x
}

to_gr_from_dt <- function(dt) GRanges(seqnames = dt$seqnames,
                                      ranges   = IRanges(start = dt$start, end = dt$end))

split_by_sample <- function(dt, type = c("amp","del")) {
  type <- match.arg(type)
  if (type == "amp") dt <- dt[copyNumber > amp_threshold] else dt <- dt[copyNumber < del_threshold]
  if (nrow(dt) == 0L) return(list())
  split(dt, by = "name", keep.by = FALSE)
}

to_ucsc <- function(chr) {
  chr <- as.character(chr)
  if (!grepl("^chr", chr[1])) chr <- paste0("chr", chr)
  chr
}

amp_threshold    <- 2                   # copyNumber > this = amplification
del_threshold    <- 2                   # copyNumber < this = deletion
bin_size         <- 1e6                   # 1 Mb bins
genome           <- "hg38"  
pdf_width        <- 8                  # inches
pdf_height       <- 8                  # inches
merge_gap        <- 1e6
std_chr <- paste0("chr", 1:22)

# Cytoband, chrom sizes, Seqinfo template 
cyt <- as.data.table(circlize::read.cytoband(species = genome))
stopifnot(ncol(cyt) >= 5)
setnames(cyt, old = names(cyt)[1:5], new = c("chr","start","end","name","gieStain"))
cyt <- cyt[chr %in% std_chr]
chrom_info <- cyt[, .(chromEnd = max(end, na.rm = TRUE)), by = .(chr)]
chrom_info <- chrom_info[.(std_chr), on = .(chr), nomatch = 0]  # preserve order

si_template <- Seqinfo(
  seqnames   = chrom_info$chr,
  seqlengths = setNames(chrom_info$chromEnd, chrom_info$chr),
  genome     = genome
)

cnvs_list <- readRDS("cnvs_list.rds") # this is generated after running MutationalTimeR.R script

cohort_dt2 <- do.call(rbind, lapply(names(cnvs_list), function(gr_name) {
  gr <- cnvs_list[[gr_name]]
  df <- as.data.table(gr)
  df$name <- gr_name  
  return(df)
}))

cohort_dt2[, seqnames := to_ucsc(seqnames)]

amps_by_sample_raw <- split_by_sample(cohort_dt2, "amp")
dels_by_sample_raw <- split_by_sample(cohort_dt2, "del")


# Rebuild BOTH lists over the same 21 samples (missing -> empty)
if (!exists("sample_ids") || is.null(sample_ids)) {
  sample_ids <- unique(cohort_dt2$name)
}

amps_by_sample <- setNames(lapply(sample_ids, function(sid)
  process_one(amps_by_sample_raw[[sid]], merge_gap, si_template)), sample_ids)

dels_by_sample <- setNames(lapply(sample_ids, function(sid)
  process_one(dels_by_sample_raw[[sid]], merge_gap, si_template)), sample_ids)


# Make fixed-size bins 

make_bins <- function(chrom_info, bin_size = 1e6) {
  gr_list <- lapply(seq_len(nrow(chrom_info)), function(i) {
    chr <- chrom_info$chr[i]
    L   <- chrom_info$chromEnd[i]
    starts <- seq.int(1, L, by = bin_size)
    ends   <- pmin(starts + bin_size - 1, L)
    GRanges(seqnames = chr, ranges = IRanges(start = starts, end = ends))
  })
  do.call(c, gr_list)
}
gr_bins <- make_bins(chrom_info, bin_size = bin_size)

# Fractions per bin (denominator = 21) 

frac_with_overlap <- function(by_sample_list, bins) {
  hits <- vapply(by_sample_list, function(gr) {
    if (length(gr) == 0L) return(rep(FALSE, length(bins)))
    countOverlaps(bins, gr) > 0
  }, logical(length(bins)))
  if (is.vector(hits)) hits <- matrix(hits, ncol = 1)
  rowMeans(hits)
}
amp_frac <- frac_with_overlap(amps_by_sample, gr_bins)
del_frac <- frac_with_overlap(dels_by_sample, gr_bins)

# Prep data frames for circos.genomicTrack
bins_df <- data.frame(
  chr   = as.character(seqnames(gr_bins)),
  start = start(gr_bins),
  end   = end(gr_bins),
  amp   = amp_frac,
  del   = del_frac
)

# Plot circos 
circos.clear()
circos.par(
  start.degree = 90,            
  clock.wise   = TRUE,
  gap.after    = rep(1, length(std_chr))
)

out_pdf          <- "cohort_cnv_frequency_circos.pdf"

pdf(out_pdf, width = pdf_width, height = pdf_height)

# Ideogram
circos.initializeWithIdeogram(species = genome, chromosome.index = std_chr, labels.cex = 1.5)

# Track 1: Amplification frequency
circos.genomicTrack(
  bins_df[, c("chr","start","end","amp")],
  ylim = c(0, 1),                      # << fixed scale
  track.height = 0.15,
  panel.fun = function(region, value, ...) {
    # (optional safety) clamp to [0,1]
    value[[1]] <- pmin(pmax(value[[1]], 0), 1)
    circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                       col = "#D6272888", border = NA, ...)
    circos.axis(major.at = c(0, .5, 1), labels = FALSE)
  }
)

# Track 2: Deletion frequency
del_df <- bins_df
del_df$del_neg <- -pmin(pmax(bins_df$del, 0), 1)  # clamp then negate

circos.genomicTrack(
  del_df[, c("chr","start","end","del_neg")],
  ylim = c(-1, 0),                    # << fixed scale
  track.height = 0.15,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                       col = "#1F77B488", border = NA, ...)
    circos.axis(major.at = c(-1, -.5, 0), labels = FALSE)
  }
)

dev.off()

# Per-chromosome-arm CNA frequencies (percent of samples)

# Rebuild cytobands cleanly and derive p/q arms
cyt_dt <- as.data.table(circlize::read.cytoband(species = genome))
stopifnot(ncol(cyt_dt) >= 5)
setnames(cyt_dt, old = names(cyt_dt)[1:5],
         new = c("chr","start","end","name","gieStain"))

# Keep autosomes; keep only true p/q bands
cyt_dt <- cyt_dt[chr %in% paste0("chr", 1:22)]
cyt_dt[, arm := fifelse(startsWith(name, "p"), "p",
                        fifelse(startsWith(name, "q"), "q", NA_character_))]
cyt_dt <- cyt_dt[!is.na(arm)]

# Collapse to contiguous p/q arm ranges
arms_dt <- cyt_dt[, .(start = min(start, na.rm = TRUE),
                      end   = max(end,   na.rm = TRUE)),
                  by = .(chr, arm)]
setorder(arms_dt, chr, arm)

# GRanges for arms
arms_gr <- GRanges(seqnames = arms_dt$chr,
                   ranges   = IRanges(start = arms_dt$start, end = arms_dt$end),
                   arm      = arms_dt$arm)

# Ensure sample_ids exists
if (!exists("sample_ids") || is.null(sample_ids)) {
  sample_ids <- unique(cohort_dt$sample_id)
}

# Fraction/percent of samples with any overlap on each arm
freq_by_arm <- function(by_sample_list, arms_gr) {
  if (length(by_sample_list) == 0L) {
    return(data.table(chr = as.character(seqnames(arms_gr)),
                      arm = mcols(arms_gr)$arm,
                      n_samples = 0L,
                      freq = 0, pct = 0))
  }
  if (is.null(names(by_sample_list))) {
    stop("by_sample_list must be a named list keyed by sample IDs")
  }
  
  overlap_mat <- sapply(by_sample_list, function(gr) {
    if (length(gr) == 0) return(rep(FALSE, length(arms_gr)))
    countOverlaps(arms_gr, gr) > 0
  })
  if (is.vector(overlap_mat)) overlap_mat <- matrix(overlap_mat, ncol = 1)
  
  frac <- rowMeans(overlap_mat)     # 0..1
  pct  <- round(100 * frac, 1)
  
  data.table(chr = as.character(seqnames(arms_gr)),
             arm = mcols(arms_gr)$arm,
             n_samples = ncol(overlap_mat),
             freq = frac,
             pct = pct)
}

# Compute amp/del per-arm frequencies
amp_arm <- freq_by_arm(amps_by_sample, arms_gr)
del_arm <- freq_by_arm(dels_by_sample, arms_gr)

setnames(amp_arm, c("freq","pct"), c("amp_freq","amp_pct"))
setnames(del_arm, c("freq","pct"), c("del_freq","del_pct"))

arm_freq <- merge(amp_arm[, .(chr, arm, n_samples, amp_freq, amp_pct)],
                  del_arm[, .(chr, arm, del_freq, del_pct)],
                  by = c("chr","arm"))
setorder(arm_freq, chr, arm)

# Add absolute counts (number of samples with amp/del on each arm)
arm_freq[, n_amp := round(amp_freq * n_samples)]
arm_freq[, n_del := round(del_freq * n_samples)]

# Quick console preview of top arms by amp or del
print(arm_freq[order(-amp_pct)][1:10])
print(arm_freq[order(-del_pct)][1:10])