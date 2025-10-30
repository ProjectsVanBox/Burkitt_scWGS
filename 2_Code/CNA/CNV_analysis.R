################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to generate circos plot
# Author: Alexander Steemers
# Date: July 2025
################################################################################

# Load required libraries

library(data.table)
library(GenomicRanges)
library(circlize)
library(biovizBase)
library(IRanges)
library(rtracklayer)
library(igraph)
library(ComplexHeatmap)

# Set output directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CNA/Figures/")

# Configurable parameters

input_dir        <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/"
genome           <- "hg38"               
bin_size         <- 1e6                   # 1 Mb bins
amp_threshold    <- 2.4                   # copyNumber > this = amplification
del_threshold    <- 1.6                   # copyNumber < this = deletion
out_pdf          <- "cohort_cnv_circos.pdf"
pdf_width        <- 8                  # inches
pdf_height       <- 8                  # inches

# Load all TSVs from purple output into a list

# Load all TSVs
tsv_files <- list.files(input_dir, pattern = "\\.integrated.cnvs.txt$", 
                        full.names = TRUE, recursive = TRUE)
# Remove unwanted files (must contain "batch", not "old") 
tsv_files_filtered <- tsv_files[grepl("batch", tsv_files) & !grepl("old", tsv_files)] 

# Load metadata 
input_df <- readxl::read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
diagnostic_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_overview.csv') 
low_callable_df<- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv') 
below_curve_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv') 
bad_baf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv') 
fail_vaf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt') 

#Make blacklist 
blacklist <- unique(c(below_curve_df$Sample_name, low_callable_df$Sample_name, bad_baf_df$Sample_name, fail_vaf_df$samplename)) 
# Extract sample IDs from filenames 
sample_ids <- sub(".integrated.cnvs.txt$", "", basename(tsv_files_filtered)) 

# Keep only non-blacklisted 
tsv_files_clean <- tsv_files_filtered[!sample_ids %in% blacklist]

samples_with_no_cnas <- c(
  "PIA9GTDBBC69", "PIA9GTDBBC75", "PIA9GTDBBC71", "PIA9GTDBBC73",
  "PIA9GTDBBC76", "PJBUGTDBBC72", "PJBUGTDABC37", "PJBUGTDBBC73",
  "PB08410-BLBM-BCELLP2C8", "PB08410-BLBM-BCELLP5G10", "PRN4GPDLBC15",
  "PVA9GTDABC42", "PVA9GTDBBC73"
)

sample_ids <- sub("\\.integrated\\.cnvs\\.txt$", "", basename(tsv_files_clean))

tsv_files_clean_2 <- tsv_files_clean[!sample_ids %in% samples_with_no_cnas]

# Function to read in tsv files

read_one <- function(f) {
  dt <- data.table::fread(f, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)
  
  # normalize headers to lowercase
  data.table::setnames(dt, tolower(names(dt)))
  
  # normalize chromosome column name
  if ("chromosome" %in% names(dt)) data.table::setnames(dt, "chromosome", "chrom")
  if ("chr" %in% names(dt) && !"chrom" %in% names(dt)) data.table::setnames(dt, "chr", "chrom")
  
  # ensure we have a CopyNumber-like column (copynumber/type/state)
  if (!"copynumber" %in% names(dt)) {
    if ("type" %in% names(dt)) {
      dt[, copynumber := type]
    } else if ("state" %in% names(dt)) {
      dt[, copynumber := state]
    } else {
      stop(sprintf("Missing CopyNumber/type/state column in %s", basename(f)))
    }
  }
  
  # required columns
  req <- c("chrom","start","end","copynumber")
  miss <- setdiff(req, names(dt))
  if (length(miss)) {
    stop(sprintf("Missing expected columns in %s: %s",
                 basename(f), paste(miss, collapse = ", ")))
  }
  
  # keep only useful columns (optional extras if present)
  keep <- intersect(c("chrom","start","end","copynumber","width","strand","rd","baf"), names(dt))
  dt <- dt[, ..keep]
  
  # rename to final casing expected downstream
  if ("copynumber" %in% names(dt)) data.table::setnames(dt, "copynumber", "CopyNumber")
  if ("rd" %in% names(dt))         data.table::setnames(dt, "rd", "RD")
  if ("baf" %in% names(dt))        data.table::setnames(dt, "baf", "BAF")
  
  # add sample id
  sid <- sub("\\.integrated\\.cnvs\\.txt$", "", basename(f), ignore.case = TRUE)
  dt[, Sample_ID := sid]
  
  dt[]
}

message(sprintf("Reading %d TSV files...", length(tsv_files_clean_2)))
cohort_dt <- rbindlist(lapply(tsv_files_clean_2, read_one), use.names = TRUE, fill = TRUE)

# Clean chromosome names to UCSC format (chr1..chr22, chrX, chrY)

to_ucsc <- function(chr) {
  chr <- as.character(chr)
  if (!grepl("^chr", chr[1])) chr <- paste0("chr", chr)
  chr
}
cohort_dt[, chrom := to_ucsc(chrom)]

# Filter to standard autosomes chromosomes present in the ideogram

std_chr <- paste0("chr", 1:22)
cohort_dt <- cohort_dt[chrom %in% std_chr & end > start]

# Ensure types
cohort_dt[, `:=`(
  start = as.numeric(start),
  end   = as.numeric(end)
)]

# merge cnvs

merge_cn_segments_iter <- function(dt,
                                   max_gap  = 30e6,
                                   min_size = 10000,
                                   n_passes = 50) {
  stopifnot(all(c("Sample_ID","chrom","start","end","CopyNumber") %in% names(dt)))
  dt <- data.table::copy(dt)
  data.table::setDT(dt)
  
  # Coerce coords; precompute width
  if (!is.integer(dt$start)) dt[, start := as.integer(start)]
  if (!is.integer(dt$end))   dt[, end   := as.integer(end)]
  dt[, width := end - start + 1L]
  
  # Optional weighted numerators/denominators for RD/BAF
  if ("RD"  %in% names(dt))  dt[, `:=`(rd_num  = RD  * width, rd_den  = width)]
  if ("BAF" %in% names(dt))  dt[, `:=`(baf_num = BAF * width, baf_den = width)]
  if (!"n_merged" %in% names(dt)) dt[, n_merged := 1L]
  
  # One merge pass within one (Sample_ID, chrom, CopyNumber)
  merge_once <- function(x, sid, chr, cn, max_gap) {
    data.table::setorder(x, start, end)
    x[, prev_end := data.table::shift(end, n = 1L, type = "lag")]
    x[, gap := start - prev_end]
    x[is.na(gap), gap := max_gap + 1L]      # force first row to start a new group
    x[, grp := cumsum(gap > max_gap)]
    
    out <- x[
      , .(
        Sample_ID  = sid,
        chrom      = chr,
        CopyNumber = cn,
        start      = min(start),
        end        = max(end),
        width      = max(end) - min(start) + 1L,
        rd_num     = if ("rd_num"  %in% names(.SD)) sum(rd_num,  na.rm = TRUE) else NA_real_,
        rd_den     = if ("rd_den"  %in% names(.SD)) sum(rd_den,  na.rm = TRUE) else NA_real_,
        baf_num    = if ("baf_num" %in% names(.SD)) sum(baf_num, na.rm = TRUE) else NA_real_,
        baf_den    = if ("baf_den" %in% names(.SD)) sum(baf_den, na.rm = TRUE) else NA_real_,
        n_merged   = sum(n_merged, na.rm = TRUE)
      )
      , by = grp
    ][, grp := NULL][]
    out
  }
  
  # ---- Iterative merging passes ----
  for (k in seq_len(n_passes)) {
    data.table::setorder(dt, Sample_ID, chrom, CopyNumber, start, end)
    dt <- dt[
      , merge_once(
        data.table::copy(.SD),
        sid = .BY$Sample_ID, chr = .BY$chrom, cn = .BY$CopyNumber, max_gap = max_gap
      )
      , by = .(Sample_ID, chrom, CopyNumber)
    ]
  }
  
  # Final RD/BAF weighted means
  if (all(c("rd_num","rd_den") %in% names(dt))) {
    dt[, RD_mean := fifelse(is.finite(rd_den) & rd_den > 0, rd_num / rd_den, NA_real_)]
  }
  if (all(c("baf_num","baf_den") %in% names(dt))) {
    dt[, BAF_mean := fifelse(is.finite(baf_den) & baf_den > 0, baf_num / baf_den, NA_real_)]
  }
  
  # Filter by size and tidy
  dt <- dt[width >= min_size]
  drop_cols <- intersect(c("rd_num","rd_den","baf_num","baf_den","prev_end","gap"), names(dt))
  if (length(drop_cols)) dt[, (drop_cols) := NULL]
  
  data.table::setorder(dt, Sample_ID, chrom, CopyNumber, start, end)
  dt[]
}

merged_dt <- merge_cn_segments_iter(
  cohort_dt,
  max_gap = 30e6, 
  min_size = 10000, 
  n_passes = 50     
)

merged_dt_filt <- merged_dt[width >= 10000]


# Split cohort_dt_merged into 6 groups by sample prefix
group1 <- merged_dt_filt[grepl("^P3G6|^PB11197", Sample_ID)]
group2 <- merged_dt_filt[grepl("^PRN4|^PB08410", Sample_ID)]
group3 <- merged_dt_filt[grepl("^P856|^PB14458", Sample_ID)]
group4 <- merged_dt_filt[grepl("^PVA9", Sample_ID)]
group5 <- merged_dt_filt[grepl("^PIA9", Sample_ID)]
group6 <- merged_dt_filt[grepl("^PJBU", Sample_ID)]


# List your groups here:
group_list <- list(group1, group2, group3, group4, group5, group6)

# Cytoband file, load once
cyto <- fread("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/cytoBandhg38.txt", 
              col.names = c("chrom", "start", "end", "band", "gieStain"))

# Function to add cytobands (same code you already tested & works)
get_cytoband_range <- function(chr, start_pos, end_pos, cyto) {
  cyto_chr <- cyto[chrom == chr]
  band_start <- cyto_chr[start <= start_pos & end >= start_pos, band]
  band_end   <- cyto_chr[start <= end_pos   & end >= end_pos, band]
  band_start <- if (length(band_start)) band_start[1L] else NA_character_
  band_end   <- if (length(band_end))   band_end[1L]   else NA_character_
  if (is.na(band_start) || is.na(band_end)) return(NA_character_)
  if (band_start == band_end) paste0(chr, band_start) else paste0(chr, band_start, "-", band_end)
}

# Store results here:
results_list <- list()

for (i in seq_along(group_list)) {
  dt <- copy(group_list[[i]])
  required_cols <- c("Sample_ID", "chrom", "CopyNumber", "start", "end", "width")
  dt <- dt[, ..required_cols]
  dt[, row_id := .I]
  dt[, `:=`(start = as.integer(start), end = as.integer(end))]
  setkey(dt, chrom, CopyNumber, start, end)
  
  # Overlaps across different samples
  ov <- foverlaps(dt, dt, type = "any", nomatch = 0L)[Sample_ID != i.Sample_ID]
  if (nrow(ov) == 0) next
  
  # Reciprocal overlap ≥ 80%
  ov[, overlap_start := pmax(start, i.start)]
  ov[, overlap_end   := pmin(end,   i.end)]
  ov[, ov_width := pmax(0L, overlap_end - overlap_start + 1L)]
  ov[, recip_x := ov_width / (end - start + 1L)]
  ov[, recip_y := ov_width / (i.end - i.start + 1L)]
  ov <- ov[recip_x >= 0.8 & recip_y >= 0.8]
  if (nrow(ov) == 0) { message("No ≥80% shared CNVs in group ", i); next }
  
  # Graph of overlaps → components
  g <- graph_from_data_frame(
    unique(ov[, .(from = row_id, to = i.row_id)]),
    directed = FALSE,
    vertices = dt[, .(row_id, Sample_ID, chrom, CopyNumber, start, end)]
  )
  comp <- components(g)
  dt[, comp_id := comp$membership[as.character(row_id)]]
  
  # --- NEW: which samples are in each component?
  comp_samples <- dt[!is.na(comp_id),
                     .(samples = list(sort(unique(Sample_ID))),   # list-column of sample IDs
                       n_samples = uniqueN(Sample_ID)),
                     by = comp_id
  ]
  
  # Keep only components shared in ≥4 samples
  keep_ids <- comp_samples[n_samples >= 2, comp_id]
  dt_filtered <- dt[comp_id %in% keep_ids]
  if (nrow(dt_filtered) == 0) next
  
  # Summarize region per component
  comp_ranges <- dt_filtered[, .(
    chrom      = unique(chrom)[1],
    CopyNumber = unique(CopyNumber)[1],
    min_start  = min(start),
    max_end    = max(end),
    width_bp   = max(end) - min(start) + 1L,
    width_mb   = round((max(end) - min(start) + 1) / 1e6, 2)
  ), by = comp_id]
  
  # Join the sample membership + counts
  comp_ranges <- comp_samples[comp_ranges, on = "comp_id"]
  
  # Cytobands
  comp_ranges[, cytoband := mapply(
    get_cytoband_range, chrom, min_start, max_end, MoreArgs = list(cyto = cyto)
  )]
  
  # Save
  results_list[[paste0("group", i)]] <- comp_ranges[]
}

# Access results
results_list$group1
results_list$group2
results_list$group3
results_list$group4
results_list$group5
results_list$group6

PVA9_samples <- c(
  "PVA9GTDBBC77","PVA9GTDBBC65","PVA9GTDBBC76","PVA9GTDBBC74","PVA9GTDBBC68","PVA9GTDBBC63",
  "PVA9GTDBBC67","PVA9GTDBBC70","PVA9GTDBBC69","PVA9GTDBBC73","PVA9GTDBBC66","PVA9GTDABC55",
  "PVA9GTDABC57","PVA9GTDABC54","PVA9GTDABC48","PVA9GTDABC21","PVA9GTDABC15","PVA9GTDABC1",
  "PVA9GTDABC52","PVA9GTDABC43","PVA9GTDABC2","PVA9GTDABC24","PVA9GTDABC34","PVA9GTDABC61",
  "PVA9GTDABC60","PVA9GTDABC59","PVA9GTDABC37","PVA9GTDABC56","PVA9GTDABC33","PVA9GTDABC53",
  "PVA9GTDABC49","PVA9GTDABC32","PVA9GTDABC40","PVA9GTDABC42","PVA9GTDABC44","PVA9GTDABC35",
  "PVA9GTDABC51","PVA9GTDABC45","PVA9GTDABC39","PVA9GTDABC50","PVA9GTDABC36","PVA9GTDABC62",
  "PVA9GTDABC58","PVA9GTDABC46","PVA9GTDABC31","PVA9GTDABC47","PVA9GTDABC18","PVA9GTDABC25",
  "PVA9GTDABC4","PVA9GTDABC16","PVA9GTDABC28","PVA9GTDABC27","PVA9GTDABC11","PVA9GTDABC19",
  "PVA9GTDABC12","PVA9GTDABC13","PVA9GTDABC23","PVA9GTDABC17","PVA9GTDABC14","PVA9GTDABC10",
  "PVA9GTDABC6","PVA9GTDABC5","PVA9GTDABC3"
)

PJBU_samples <- c(
  "PJBUGTDBBC69","PJBUGTDBBC74","PJBUGTDBBC67","PJBUGTDBBC64","PJBUGTDBBC72",
  "PJBUGTDBBC68","PJBUGTDBBC73","PJBUGTDBBC59","PJBUGTDBBC65","PJBUGTDBBC76",
  "PJBUGTDBBC71","PJBUGTDBBC81","PJBUGTDBBC70","PJBUGTDBBC63","PJBUGTDBBC80",
  "PJBUGTDBBC79","PJBUGTDABC9","PJBUGTDABC43","PJBUGTDABC54","PJBUGTDABC26",
  "PJBUGTDABC37","PJBUGTDABC24","PJBUGTDABC11","PJBUGTDABC8","PJBUGTDABC6",
  "PJBUGTDABC19","PJBUGTDABC42","PJBUGTDABC50","PJBUGTDABC13","PJBUGTDABC29",
  "PJBUGTDABC30","PJBUGTDABC33","PJBUGTDABC32","PJBUGTDABC51","PJBUGTDABC56",
  "PJBUGTDABC3","PJBUGTDABC40","PJBUGTDABC27","PJBUGTDABC55","PJBUGTDABC16",
  "PJBUGTDABC2","PJBUGTDABC31","PJBUGTDABC44","PJBUGTDABC38","PJBUGTDABC39",
  "PJBUGTDABC21","PJBUGTDABC58","PJBUGTDABC53","PJBUGTDABC49","PJBUGTDABC52",
  "PJBUGTDABC1","PJBUGTDABC18","PJBUGTDABC45","PJBUGTDABC17","PJBUGTDABC41",
  "PJBUGTDABC23","PJBUGTDABC35","PJBUGTDABC25"
)

P856_samples <-  c( "P856GDDUBC42", "P856GDDUBC44", "PB14458-BLPL-BCELLP4L3", "PB14458-BLPL-BCELLP4K5", "PB14458-BLPL-BCELLP4M3", "PB14458-BLPL-BCELLP4B3", "PB14458-BLPL-BCELLP4L5" , "P856GDDBBC63", "PB14458-BLBM-BCELLP2F2", "PB14458-BLBM-BCELLP2B3", "P856GDDBBC60", "P856GDDBBC46", "PB14458-BLBM-BCELLP2B4", "PB14458-BLBM-BCELLP2C4", "PB14458-BLBM-BCELLP2N2", "P856GDDBBC59", "P856GDDBBC64", "PB14458-BLBM-BCELLP2F4", "P856GDDBBC58", "P856GDDBBC54", "P856GDDBBC48", "PB14458-BLBM-BCELLP2I2", "PB14458-BLBM-BCELLP2L4",  "PB14458-BLBM-BCELLP2E4", "P856GDDBBC57", "P856GDDBBC62" ,"PB14458-BLBM-BCELLP2L3", "P856GDDBBC61") 

P3G6_samples <- c(
  "P3G6GPDABC31", "PB11197-BLASC-BCELLP2F4", "PB11197-BLASC-BCELLP2D4", 
  "PB11197-BLASC-BCELLP2C4", "PB11197-BLASC-BCELLP2B4", "PB11197-BLASC-BCELLP2E4", 
  "P3G6GPDABC28", "PB11197-BLASC-BCELLP1P3", "PB11197-BLASC-BCELLP1C4", 
  "PB11197-BLASC-BCELLP1L3", "PB11197-BLASC-BCELLP1J3", "PB11197-BLASC-BCELLP1K4", 
  "PB11197-BLASC-BCELLP1O3", "PB11197-BLASC-BCELLP1I4", "PB11197-BLASC-BCELLP1B4", 
  "P3G6GPDABC26"
)

PIA9_samples <- c(
  "PIA9GTDBBC75","PIA9GTDBBC64","PIA9GTDBBC52","PIA9GTDBBC54","PIA9GTDBBC57",
  "PIA9GTDBBC73","PIA9GTDBBC58","PIA9GTDBBC60","PIA9GTDBBC67","PIA9GTDBBC59",
  "PIA9GTDBBC63","PIA9GTDBBC72","PIA9GTDBBC61","PIA9GTDBBC77","PIA9GTDBBC68",
  "PIA9GTDBBC55","PIA9GTDBBC65","PIA9GTDBBC53","PIA9GTDBBC66","PIA9GTDBBC56",
  "PIA9GTDBBC74","PIA9GTDBBC71","PIA9GTDBBC76",
  "PIA9GTDABC37","PIA9GTDABC7","PIA9GTDABC49","PIA9GTDABC21","PIA9GTDABC33",
  "PIA9GTDABC19","PIA9GTDABC9","PIA9GTDABC26","PIA9GTDABC44","PIA9GTDABC5",
  "PIA9GTDABC6","PIA9GTDABC25","PIA9GTDABC48","PIA9GTDABC18","PIA9GTDABC40",
  "PIA9GTDABC34","PIA9GTDABC50","PIA9GTDABC42","PIA9GTDABC46","PIA9GTDABC47",
  "PIA9GTDABC43","PIA9GTDABC20","PIA9GTDABC36","PIA9GTDABC24","PIA9GTDABC39",
  "PIA9GTDABC15","PIA9GTDABC51","PIA9GTDABC23","PIA9GTDABC17","PIA9GTDABC45",
  "PIA9GTDABC41","PIA9GTDABC27","PIA9GTDABC35","PIA9GTDABC22"
)

PRN4_samples <- c("PB08410-BLBM-BCELLP5G10", "PRN4GPDBBC07", "PB08410-BLBM-BCELLP2C8", "PB08410-BLBM-BCELLP2G8", "PB08410-BLBM-BCELLP2D8", "PB08410-BLBM-BCELLP5E8", "PB08410-BLBM-BCELLP5F8", "PRN4GPDLBC15", "PRN4GPDLBC17","PRN4GPDLBC09","PRN4GPDLBC11", "PB08410-BLLN-BCELLP4F10",  "PRN4GPDLBC21", "PRN4GPDLBC23", "PB08410-BLLN-BCELLP4D10", "PB08410-BLLN-BCELLP2D10","PRN4GPDLBC22", "PB08410-BLLN-BCELLP4B11", "PB08410-BLLN-BCELLP1B11", "PRN4GPDLBC20", "PB08410-BLLN-BCELLP2B10","PB08410-BLLN-BCELLP2E10",  "PRN4GPDLBC10", "PB08410-BLLN-BCELLP4G10", "PRN4GPDLBC16", "PRN4GPDLBC19")


# ----- your cohort vectors defined above (P3G6_samples, P856_samples, ... PJBU_samples) -----

sample_sets <- list(
  P3G6 = P3G6_samples,
  P856 = P856_samples,
  PIA9 = PIA9_samples,
  PRN4 = PRN4_samples,
  PVA9 = PVA9_samples,
  PJBU = PJBU_samples
)

# ---- Step 1: Merge all groups (unchanged) ----
all_cnvs <- rbindlist(
  lapply(names(results_list), function(nm) {
    df <- results_list[[nm]]
    if (is.null(df) || !nrow(df)) return(NULL)
    df[, group := nm]
    df
  }),
  use.names = TRUE, fill = TRUE
)
stopifnot(!is.null(all_cnvs), "No CNVs found in results_list")

# ---- Step 2: Build SIGNED matrix with LOH (=2) ----
# Expand samples; keep CopyNumber so we can map to categories
long_signed <- rbindlist(lapply(results_list, function(df) {
  if (is.null(df) || !nrow(df)) return(NULL)
  df[, .(cytoband, CopyNumber, Sample_ID = unlist(samples)), by = comp_id]
}), fill = TRUE)

# Map CopyNumber to codes: +1 Gain, -1 Loss, 2 LOH, 0 other/unknown
long_signed[, code :=
              fifelse(grepl("(?i)\\b(gain|amp)\\b", CopyNumber),  1L,
                      fifelse(grepl("(?i)\\b(loss|del)\\b", CopyNumber), -1L,
                              fifelse(grepl("(?i)\\b(loh|cnn[-_ ]?loh|upd)\\b", CopyNumber), 2L, 0L)))
]

# Collapse duplicates per (cytoband, Sample_ID)
# Priority rule: Gain (1) > Loss (-1) > LOH (2) > None (0)
# Implement by checking presence in that order
long_signed <- long_signed[
  , .(code = if (any(code == 1L)) 1L
      else if (any(code == -1L)) -1L
      else if (any(code == 2L)) 2L
      else 0L),
  by = .(cytoband, Sample_ID)
]

# Wide matrix
wide_signed <- dcast(long_signed, cytoband ~ Sample_ID,
                     value.var = "code", fun.aggregate = max, fill = 0)
signed_mat_all <- as.matrix(wide_signed[, -1, with = FALSE])
rownames(signed_mat_all) <- wide_signed$cytoband

# ---- Helper: conform columns to a given order (add missing as zeros) ----
conform_to_subset <- function(mat, wanted_cols) {
  keep <- intersect(wanted_cols, colnames(mat))
  mat_sub <- if (length(keep)) mat[, keep, drop = FALSE] else
    matrix(0L, nrow = nrow(mat), ncol = 0, dimnames = list(rownames(mat), NULL))
  missing <- setdiff(wanted_cols, colnames(mat_sub))
  if (length(missing)) {
    add <- matrix(0L, nrow = nrow(mat_sub), ncol = length(missing),
                  dimnames = list(rownames(mat_sub), missing))
    mat_sub <- if (ncol(mat_sub) == 0L) add else cbind(mat_sub, add)
  }
  mat_sub[, wanted_cols, drop = FALSE]
}

# ---- Step 3: Plot one heatmap per cohort ----
# Discrete palette & breaks for {-1, 0, 1, 2} = {Loss, None, Gain, LOH}
heat_colors <- c("blue", "white", "red", "green")
# Breaks length must be colors length + 1; set bins centered on the integers
heat_breaks <- c(-1.5, -0.5, 0.5, 1.5, 2.5)

for (set_name in names(sample_sets)) {
  wanted <- sample_sets[[set_name]]
  # subset columns to cohort (add blanks for missing)
  mat_set <- conform_to_subset(signed_mat_all, wanted)
  
  # keep only rows present in this cohort (any non-zero, including LOH=2)
  present_rows <- rowSums(mat_set != 0) > 0
  mat_set <- mat_set[present_rows, , drop = FALSE]
  if (nrow(mat_set) == 0L) {
    message("No cytobands present for cohort ", set_name, " — skipping.")
    next
  }
  
  # dynamic height so row names don't crush
  pdf_height <- max(6, min(20, nrow(mat_set) * 0.20))
  
  pheatmap::pheatmap(
    mat_set,
    cluster_rows = TRUE,
    cluster_cols = FALSE,              # keep your cohort order
    color = heat_colors,
    breaks = heat_breaks,
    main = paste0("CNV Gains/Losses/LOH – ", set_name, " (all groups combined)"),
    fontsize = 8,
    show_colnames = TRUE,
    fontsize_col = 6,
    angle_col = 90,
    legend_labels = c("Loss", "None", "Gain", "LOH"),
    filename = paste0("CNV_heatmap_", set_name, "_gains_losses_loh.pdf"),
    width = 10, height = pdf_height
  )
  
  # Save matrix used
  out_df <- cbind(cytoband = rownames(mat_set), as.data.frame(mat_set, check.names = FALSE))
  write.csv(out_df, paste0("CNV_matrix_", set_name, "_gains_losses_loh_present_only.csv"), row.names = FALSE)
}

###########################################################


# # Keep only CNAs ≥ 10 Mb
# cohort_dt_filt <- cohort_dt[width >= 10e6]
# 
# # Per-sample GRanges for amp & del
# 
# make_gr <- function(dt, type = c("amp","del","loh")) {
#   type <- match.arg(type)
#   # map type -> label in CopyNumber
#   lab <- switch(type,
#                 amp = "Gain",
#                 del = "Loss",
#                 loh = "LOH")
#   
#   # filter by the label
#   dt <- dt[CopyNumber == lab]
#   
#   if (nrow(dt) == 0) return(list())
#   
#   # split by the right sample column
#   split(dt, by = "Sample_ID", keep.by = FALSE)
# }
# 
# amps_by_sample <- make_gr(cohort_dt, "amp")
# dels_by_sample <- make_gr(cohort_dt, "del")
# 
# # Merge CNAs of the same type if they're <10 Mb apart
# 
# merge_gap <- 10e6
# 
# to_gr <- function(x) GRanges(seqnames = x$chrom, ranges = IRanges(start = x$start, end = x$end))
# 
# amps_by_sample <- lapply(amps_by_sample, function(dt) {
#   if (nrow(dt) == 0) return(GRanges())
#   reduce(to_gr(dt), min.gapwidth = merge_gap)
# })
# 
# dels_by_sample <- lapply(dels_by_sample, function(dt) {
#   if (nrow(dt) == 0) return(GRanges())
#   reduce(to_gr(dt), min.gapwidth = merge_gap)
# })
# 
# sample_ids <- unique(cohort_dt$sample_id)
# N <- length(sample_ids)
# 
# # Read cytoband
# 
# cyt <- circlize::read.cytoband(species = genome)
# cyt <- as.data.table(cyt)
# 
# # Normalize cyt columns
# 
# setnames(cyt, old = c("df.V1","df.V2","df.V3"), new = c("chr","start","end"), skip_absent = TRUE)
# 
# # Per-chromosome lengths
# 
# chrom_info <- cyt[, .(chromEnd = max(end, na.rm = TRUE)), by = .(chr)]
# 
# # Keep/order standard chromosomes (this preserves std_chr order)
# 
# std_chr <- paste0("chr", 1:22)
# chrom_info <- chrom_info[.(std_chr), on = .(chr), nomatch = 0]
# 
# # (optional) sanity check
# stopifnot(identical(chrom_info$chr, std_chr[std_chr %in% chrom_info$chr]))
# 
# # Make fixed-size bins
# 
# make_bins <- function(chrom_info, bin_size = 1e6) {
#   gr_list <- lapply(seq_len(nrow(chrom_info)), function(i) {
#     chr <- chrom_info$chr[i]
#     L   <- chrom_info$chromEnd[i]
#     starts <- seq.int(1, L, by = bin_size)
#     ends   <- pmin(starts + bin_size - 1, L)
#     GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = starts, end = ends))
#   })
#   do.call(c, gr_list)
# }
# 
# gr_bins <- make_bins(chrom_info, bin_size = bin_size)
# 
# # Frequency calculation (fraction of samples with any overlap per bin)
# 
# frac_with_overlap <- function(by_sample_list, bins) {
#   if (length(by_sample_list) == 0) return(rep(0, length(bins)))
#   hits_per_sample <- vapply(by_sample_list, function(gr) {
#     if (length(gr) == 0) return(rep(FALSE, length(bins)))
#     countOverlaps(bins, gr) > 0
#   }, logical(length(bins)))
#   # if only 1 sample had events, vapply returns vector; coerce to matrix
#   if (is.vector(hits_per_sample)) hits_per_sample <- matrix(hits_per_sample, ncol = 1)
#   rowMeans(hits_per_sample) # fraction of samples
# }
# 
# amp_frac <- frac_with_overlap(amps_by_sample, gr_bins)
# del_frac <- frac_with_overlap(dels_by_sample, gr_bins)
# 
# # Prepare per-chromosome data frames for plotting with circlize
# 
# bins_df <- data.frame(
#   chr   = as.character(seqnames(gr_bins)),
#   start = start(gr_bins),
#   end   = end(gr_bins),
#   amp   = amp_frac,
#   del   = del_frac
# )
# 
# # Plotting
# 
# circos.clear()
# circos.par(
#   start.degree = 90,    # chr1 starts at 12 o'clock
#   clock.wise   = TRUE,  # proceed clockwise
#   gap.after    = rep(1, length(std_chr))  # uniform small gaps
# )
# pdf(out_pdf, width = pdf_width, height = pdf_height)
# circos.initializeWithIdeogram(species = genome, chromosome.index = std_chr, labels.cex = 1.5)
# 
# # Track 1: Amplification frequency (outer ring)
# 
# max_y <- max(0.01, max(bins_df$amp, na.rm = TRUE))
# circos.genomicTrack(
#   bins_df[, c("chr","start","end","amp")],
#   ylim = c(0, max_y),
#   track.height = 0.15,
#   panel.fun = function(region, value, ...) {
#     circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
#                        col = "#D6272888", border = NA, ...)  # red-ish, semi-transparent
#     circos.axis(major.at = NULL, labels = FALSE)
#   }
# )
# 
# # Track 2: Deletion frequency (inner ring; plotted negative so it sits inside)
# 
# min_y <- -max(0.01, max(bins_df$del, na.rm = TRUE))
# del_df <- bins_df
# del_df$del_neg <- -del_df$del
# circos.genomicTrack(
#   del_df[, c("chr","start","end","del_neg")],
#   ylim = c(min_y, 0),
#   track.height = 0.15,
#   panel.fun = function(region, value, ...) {
#     circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
#                        col = "#1F77B488", border = NA, ...)  # blue-ish, semi-transparent
#     circos.axis(major.at = NULL, labels = FALSE)
#   }
# )
# 
# dev.off()
# 
# 
# # Helper: build circos-ready segments for a single sample
# .build_seg_df <- function(dt_one) {
#   if (nrow(dt_one) == 0) return(NULL)
#   dt_one <- dt_one[end > start & chrom %in% paste0("chr", 1:22)]
#   if (nrow(dt_one) == 0) return(NULL)
#   
#   # Map CopyNumber labels to colors (Gain / Loss / LOH)
#   col_map <- c(Gain = "#E64B35", Loss = "#4C78A8", LOH = "#56B870")
#   dt_one <- dt_one[CopyNumber %in% names(col_map)]
#   if (nrow(dt_one) == 0) return(NULL)
#   
#   data.frame(
#     chr   = dt_one$chrom,
#     start = as.numeric(dt_one$start),
#     end   = as.numeric(dt_one$end),
#     type  = as.character(dt_one$CopyNumber),
#     stringsAsFactors = FALSE
#   )
# }
# 
# # Function: plot one sample to its own PDF
# plot_sample_circos <- function(sample_id, dt_sample,
#                                genome = "hg38",
#                                std_chr = paste0("chr", 1:22),
#                                width = 7, height = 7) {
#   seg_df <- .build_seg_df(dt_sample)
#   if (is.null(seg_df) || nrow(seg_df) == 0) {
#     message("No plottable CNAs for ", sample_id, " — skipping.")
#     return(invisible(NULL))
#   }
#   
#   # start a new PDF for this sample
#   pdf(file.path(getwd(), paste0(sample_id, "_cna_circos.pdf")),
#       width = width, height = height)
#   
#   # circos setup
#   circos.clear()
#   circos.par(start.degree = 90, gap.after = rep(1, length(std_chr)))
#   circos.initializeWithIdeogram(species = genome,
#                                 chromosome.index = std_chr,
#                                 labels.cex = 0.6)
#   
#   # one track for all CNAs, colored by type
#   circos.genomicTrack(
#     seg_df[, c("chr", "start", "end")],
#     ylim = c(0, 1), track.height = 0.12,
#     panel.fun = function(region, value, ...) {}
#   )
#   
#   # draw rectangles per type for clarity
#   col_map <- c(Gain = "#E64B35", Loss = "#4C78A8", LOH = "#56B870")
#   for (typ in names(col_map)) {
#     sub <- seg_df[seg_df$type == typ, c("chr","start","end"), drop = FALSE]
#     if (nrow(sub) == 0) next
#     circos.genomicRect(
#       sub, ytop = 1, ybottom = 0,
#       col = grDevices::adjustcolor(col_map[typ], alpha.f = 0.75),
#       border = NA
#     )
#   }
#   
#   title(paste("CNAs —", sample_id))
#   dev.off()
# }
# 
# # ---- Run for each sample ----
# samples <- sort(unique(cohort_dt$Sample_ID))
# invisible(lapply(samples, function(sid) {
#   plot_sample_circos(sid, cohort_dt[Sample_ID == sid])
# }))
# 
# 
# get_cytobands <- function(genome = c("hg38","hg19"), chrs = NULL) {
#   genome <- match.arg(genome)
#   
#   # ---- Try AnnotationHub first (stable + cached) ----
#   ah_ok <- TRUE
#   gr <- tryCatch({
#     ah <- AnnotationHub::AnnotationHub()
#     q <- AnnotationHub::query(ah, c("UCSC","cytoBand", genome))
#     # pick a GRanges resource with correct genome
#     idx <- which(vapply(q, function(x) inherits(ah[[names(x)]], "GRanges"), logical(1)))
#     if (length(idx) > 0) {
#       cand_id <- names(q)[idx[1]]
#       ah[[cand_id]]
#     } else {
#       stop("No GRanges cytoband in AnnotationHub")
#     }
#   }, error = function(e) { ah_ok <<- FALSE; NULL })
#   
#   # ---- Fallback: direct UCSC by TABLE (no 'track') ----
#   if (!ah_ok) {
#     sess <- browserSession("UCSC")
#     genome(sess) <- genome
#     for (tname in c("cytoBandIdeo","cytoBand")) {
#       tbl <- tryCatch({
#         q <- ucscTableQuery(sess, table = tname)  # <— no 'track' argument
#         getTable(q)
#       }, error = function(e) NULL)
#       if (!is.null(tbl) && all(c("chrom","chromStart","chromEnd","name") %in% names(tbl))) {
#         gr <- GRanges(
#           seqnames = tbl$chrom,
#           ranges   = IRanges(start = tbl$chromStart + 1, end = tbl$chromEnd),
#           name     = tbl$name,
#           gieStain = tbl$gieStain
#         )
#         break
#       }
#     }
#     if (is.null(gr)) stop("Couldn't retrieve cytobands from UCSC for ", genome)
#   }
#   
#   # keep UCSC-style seqnames (chr1 etc.)
#   if (!is.null(chrs)) gr <- gr[seqnames(gr) %in% chrs]
#   gr[order(as.character(seqnames(gr)), start(gr))]
# }
# 
# annotate_with_cytobands <- function(df, genome = c("hg38","hg19"),
#                                     chr_col = "chrom", start_col = "start", end_col = "end") {
#   genome <- match.arg(genome)
#   stopifnot(all(c(chr_col, start_col, end_col) %in% names(df)))
#   
#   # Use helper we already defined
#   ideo <- get_cytobands(genome, chrs = unique(df[[chr_col]]))
#   
#   qry  <- GRanges(df[[chr_col]], IRanges(df[[start_col]], df[[end_col]]))
#   hits <- findOverlaps(qry, ideo, type = "any")
#   bands_by_q <- split(mcols(ideo)$name[subjectHits(hits)], queryHits(hits))
#   
#   compact_label <- function(chr, bands) {
#     if (length(bands) == 0) return(NA_character_)
#     sb <- bands[1]; eb <- bands[length(bands)]
#     if (substr(sb,1,1) == substr(eb,1,1)) paste0(chr, sb, "–", sub("^[pq]", "", eb))
#     else paste0(chr, sb, "–", eb)
#   }
#   
#   annot <- lapply(seq_along(qry), function(i) {
#     bands <- bands_by_q[[as.character(i)]]
#     if (is.null(bands) || length(bands) == 0) {
#       data.frame(start_band = NA_character_, end_band = NA_character_,
#                  bands = NA_character_, cytoband_range = NA_character_)
#     } else {
#       data.frame(start_band = bands[1],
#                  end_band   = bands[length(bands)],
#                  bands = paste(bands, collapse = ","),
#                  cytoband_range = compact_label(as.character(seqnames(qry)[i]), bands))
#     }
#   })
#   annot_df <- do.call(rbind, annot)
#   
#   # Bind annotations to the original dataframe
#   cbind(df, annot_df)
# }
# 
# group1_with_cytoband <- annotate_with_cytobands(
#   group1,
#   genome = "hg38",
#   chr_col = "chrom",
#   start_col = "start",
#   end_col = "end"
# )
# 
# group2_with_cytoband <- annotate_with_cytobands(
#   group2,
#   genome = "hg38",
#   chr_col = "chrom",
#   start_col = "start",
#   end_col = "end"
# )
# 
# group3_with_cytoband <- annotate_with_cytobands(
#   group3,
#   genome = "hg38",
#   chr_col = "chrom",
#   start_col = "start",
#   end_col = "end"
# )
# 
# group4_with_cytoband <- annotate_with_cytobands(
#   group4,
#   genome = "hg38",
#   chr_col = "chrom",
#   start_col = "start",
#   end_col = "end"
# )
# 
# group5_with_cytoband <- annotate_with_cytobands(
#   group5,
#   genome = "hg38",
#   chr_col = "chrom",
#   start_col = "start",
#   end_col = "end"
# )
# 
# group6_with_cytoband <- annotate_with_cytobands(
#   group6,
#   genome = "hg38",
#   chr_col = "chrom",
#   start_col = "start",
#   end_col = "end"
# )
# 
# 
# 
