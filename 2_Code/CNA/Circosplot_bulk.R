################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to generate circos plot
# Author: Alexander Steemers
# Date: July 2025
################################################################################

# Load required libraries

packages <- c("data.table", "GenomicRanges", "circlize")
for (p in packages) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table)
library(GenomicRanges)
library(circlize)

# Set output directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CNA/Figures/")

# Configurable parameters

input_dir        <- "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input"
genome           <- "hg38"               
bin_size         <- 1e6                   # 1 Mb bins
amp_threshold    <- 2.4                   # copyNumber > this = amplification
del_threshold    <- 1.6                   # copyNumber < this = deletion
out_pdf          <- "cohort_cnv_frequency_circos.pdf"
pdf_width        <- 8                  # inches
pdf_height       <- 8                  # inches

# Load all TSVs from purple output into a list

tsv_files <- list.files(input_dir, pattern = ".bwamem2.samtools.gatk4spark.dedup.purple.cnv.somatic.tsv$", full.names = TRUE, recursive = TRUE)
tsv_files <- tsv_files[!grepl("old|work|wrong", tsv_files)]
tsv_files <- tsv_files[!grepl("IPY|KXA", tsv_files)] # first is copy of PVA9 and the other doesn't have MYC translocation
length(tsv_files) # 21 samples now (as expected)


# Function to read in tsv files

read_one <- function(f) {
  dt <- fread(f, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
  # Normalize column names
  setnames(dt, tolower(names(dt)))
  # Keep essentials; rename for consistency
  keep <- c("chromosome","start","end","copynumber")
  missing <- setdiff(keep, names(dt))
  if (length(missing)) stop(sprintf("Missing expected columns in %s: %s", basename(f), paste(missing, collapse = ", ")))
  dt[, sample_id := tools::file_path_sans_ext(basename(f))]
  dt[]
}

message(sprintf("Reading %d TSV files...", length(tsv_files)))
cohort_dt <- rbindlist(lapply(tsv_files, read_one), use.names = TRUE, fill = TRUE)

# Clean chromosome names to UCSC format (chr1..chr22, chrX, chrY)

to_ucsc <- function(chr) {
  chr <- as.character(chr)
  if (!grepl("^chr", chr[1])) chr <- paste0("chr", chr)
  chr
}
cohort_dt[, chromosome := to_ucsc(chromosome)]

# Filter to standard autosomes chromosomes present in the ideogram

std_chr <- paste0("chr", 1:22)
cohort_dt <- cohort_dt[chromosome %in% std_chr & end > start]

# Filter to standard autosomes + only long enough CNAs (>10 Mb)

std_chr <- paste0("chr", 1:22)   # autosomes only
cohort_dt <- cohort_dt[
  chromosome %in% std_chr & end > start & (end - start + 1) >= 10e6
]

# Per-sample GRanges for amp & del

make_gr <- function(dt, type = c("amp","del")) {
  type <- match.arg(type)
  if (type == "amp") {
    dt <- dt[copynumber > amp_threshold]
  } else {
    dt <- dt[copynumber < del_threshold]
  }
  if (nrow(dt) == 0) return(list())
  split(dt, by = "sample_id", keep.by = FALSE)
}

amps_by_sample <- make_gr(cohort_dt, "amp")
dels_by_sample <- make_gr(cohort_dt, "del")

# Merge CNAs of the same type if they're <10 Mb apart

merge_gap <- 10e6

to_gr <- function(x) GRanges(seqnames = x$chromosome, ranges = IRanges(start = x$start, end = x$end))

amps_by_sample <- lapply(amps_by_sample, function(dt) {
  if (nrow(dt) == 0) return(GRanges())
  reduce(to_gr(dt), min.gapwidth = merge_gap)
})

dels_by_sample <- lapply(dels_by_sample, function(dt) {
  if (nrow(dt) == 0) return(GRanges())
  reduce(to_gr(dt), min.gapwidth = merge_gap)
})

sample_ids <- unique(cohort_dt$sample_id)
N <- length(sample_ids)

# Read cytoband

cyt <- circlize::read.cytoband(species = genome)
cyt <- as.data.table(cyt)

# Normalize cyt columns

setnames(cyt, old = c("df.V1","df.V2","df.V3"), new = c("chr","start","end"), skip_absent = TRUE)

# Per-chromosome lengths

chrom_info <- cyt[, .(chromEnd = max(end, na.rm = TRUE)), by = .(chr)]

# Keep/order standard chromosomes (this preserves std_chr order)

std_chr <- paste0("chr", 1:22)
chrom_info <- chrom_info[.(std_chr), on = .(chr), nomatch = 0]

# (optional) sanity check
stopifnot(identical(chrom_info$chr, std_chr[std_chr %in% chrom_info$chr]))

# Make fixed-size bins

make_bins <- function(chrom_info, bin_size = 1e6) {
  gr_list <- lapply(seq_len(nrow(chrom_info)), function(i) {
    chr <- chrom_info$chr[i]
    L   <- chrom_info$chromEnd[i]
    starts <- seq.int(1, L, by = bin_size)
    ends   <- pmin(starts + bin_size - 1, L)
    GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = starts, end = ends))
  })
  do.call(c, gr_list)
}

gr_bins <- make_bins(chrom_info, bin_size = bin_size)

# Frequency calculation (fraction of samples with any overlap per bin)

frac_with_overlap <- function(by_sample_list, bins) {
  if (length(by_sample_list) == 0) return(rep(0, length(bins)))
  hits_per_sample <- vapply(by_sample_list, function(gr) {
    if (length(gr) == 0) return(rep(FALSE, length(bins)))
    countOverlaps(bins, gr) > 0
  }, logical(length(bins)))
  # if only 1 sample had events, vapply returns vector; coerce to matrix
  if (is.vector(hits_per_sample)) hits_per_sample <- matrix(hits_per_sample, ncol = 1)
  rowMeans(hits_per_sample) # fraction of samples
}

amp_frac <- frac_with_overlap(amps_by_sample, gr_bins)
del_frac <- frac_with_overlap(dels_by_sample, gr_bins)

# Prepare per-chromosome data frames for plotting with circlize

bins_df <- data.frame(
  chr   = as.character(seqnames(gr_bins)),
  start = start(gr_bins),
  end   = end(gr_bins),
  amp   = amp_frac,
  del   = del_frac
)

# Plotting

circos.clear()
circos.par(
  start.degree = 90,    # chr1 starts at 12 o'clock
  clock.wise   = TRUE,  # proceed clockwise
  gap.after    = rep(1, length(std_chr))  # uniform small gaps
)
pdf(out_pdf, width = pdf_width, height = pdf_height)
circos.initializeWithIdeogram(species = genome, chromosome.index = std_chr, labels.cex = 1.5)

# Track 1: Amplification frequency (outer ring)

max_y <- max(0.01, max(bins_df$amp, na.rm = TRUE))
circos.genomicTrack(
  bins_df[, c("chr","start","end","amp")],
  ylim = c(0, max_y),
  track.height = 0.15,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                       col = "#D6272888", border = NA, ...)  # red-ish, semi-transparent
    circos.axis(major.at = NULL, labels = FALSE)
  }
)

# Track 2: Deletion frequency (inner ring; plotted negative so it sits inside)

min_y <- -max(0.01, max(bins_df$del, na.rm = TRUE))
del_df <- bins_df
del_df$del_neg <- -del_df$del
circos.genomicTrack(
  del_df[, c("chr","start","end","del_neg")],
  ylim = c(min_y, 0),
  track.height = 0.15,
  panel.fun = function(region, value, ...) {
    circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                       col = "#1F77B488", border = NA, ...)  # blue-ish, semi-transparent
    circos.axis(major.at = NULL, labels = FALSE)
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

