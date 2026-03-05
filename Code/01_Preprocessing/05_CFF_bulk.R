################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Script to obtain clonal bulk mutations by using the cancer cell fraction method
# Author: Alexander Steemers
################################################################################

library(ccube)
library(MutationTimeR)
library(ggplot2)
library(stringr)
library(dplyr)
library(MutationalPatterns)
library(reshape2)
library(ggbeeswarm)
library(cowplot)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.NCBI.GRCh38)
ref_genome <- "BSgenome.Hsapiens.NCBI.GRCh38"

# Get smurf filtered vcf files for the bulk samples (these should be downloaded from Mendeley DOI: 10.17632/phk5jhhm7d.1)
# Find out which VFCs to download by checking the file names below
# They should be 21 bulk samples in total

bulk_file_names <- readLines("Data/Bulk_smurf_filtered_file_names.txt")

#all_filtered_vcfs <- "~/path/to/vcf/files.vcf"

samples_names <- c()
full_names <- c()

for (vcf_file in vcf_files) {
  sample_name <- sub(".*?/SMuRF/(.*)\\.vep.*", "\\1", vcf_file)
  
  full_names <- c(full_names, sample_name)
  samples_names <- c(samples_names, sample_name)
}

cnv_samps <- samples_names

# Get cnv files for the bulk samples (these should be downloaded from Mendeley DOI: 10.17632/phk5jhhm7d.1)
# Find out which files to download by checking the file names below
# They should be 21 bulk samples in total

cnv_files <- readLines("Data/cnv_files_bulk.txt")

# Get purity files for the bulk samples (these should be downloaded from Mendeley DOI: 10.17632/phk5jhhm7d.1)
# Find out which files to download by checking the file names below
# They should be 21 bulk samples in total

purity_f <- readLines("Data/purity_files_bulk.txt")

# load data
vcfs = lapply(vcf_files, function(v) {
  cat(v, fill = T)
  VariantAnnotation::readVcf(v)
})
lengths(vcfs)
cnvs = lapply(cnv_files, read.table, sep = "\t", header = T)
purities = lapply(purity_f, read.table, sep = "\t", header = T)
names(vcfs) = names(cnvs) = names(purities) = cnv_samps

tumour_ids <- c(
  "PMABM000GTJ","PMABM000DRF","PMABM000GUL","PRN4GBDLBC72","PMABM000FZO",
  "PMLBM000AUV","PB14458-BLBM-BCELLBULK","PB14458-BLPL-BCELLBULK",
  "PMABM000HDI","PMABM000HDS","PVA9GBDABC78","PMLBM000IDI","PMLBM000KOD",
  "PMABM000GLZ","PMLBM000ILF","PMABM000GPM","PMLBM000LLM","PMLBM000CVT",
  "PJBUGBDABC82","PB11197-BLASC-BCELLBULK","PIA9GBDABC78"
)
tumour_ids <- unique(tumour_ids)

vcf_to_ssm_counts_one <- function(vcf_obj, tumour_ids) {
  
  samp <- VariantAnnotation::samples(header(vcf_obj))
  tum_hits <- which(samp %in% tumour_ids)
  
  if (length(tum_hits) == 0) {
    stop("No tumour sample found in VCF. Samples are: ", paste(samp, collapse = ", "))
  }
  if (length(tum_hits) > 1) {
    stop("More than one tumour ID matched in this VCF: ", paste(samp[tum_hits], collapse = ", "))
  }
  
  tum_idx <- tum_hits[1]
  
  AD <- geno(vcf_obj)$AD
  AD_tum <- AD[, tum_idx]
  
  ref_counts <- vapply(AD_tum, function(x) as.integer(x[1]), integer(1))
  var_counts <- vapply(AD_tum, function(x) as.integer(x[2]), integer(1))
  
  rr <- rowRanges(vcf_obj)
  chr <- as.character(seqnames(rr))
  pos <- start(rr)
  
  ssm <- data.frame(
    mutation_id   = paste0(chr, "_", pos),
    chr           = chr,
    pos           = pos,
    ref_counts    = ref_counts,
    var_counts    = var_counts,
    tumour_sample = samp[tum_idx],
    stringsAsFactors = FALSE
  )
  
  subset(ssm, (ref_counts + var_counts) > 0)
}
ssm_list <- lapply(vcfs, vcf_to_ssm_counts_one, tumour_ids = tumour_ids)

get_purity_value <- function(p) {
  # p can be file path OR data.frame
  pur <- if (is.character(p) && length(p) == 1) {
    read.delim(p, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    p
  }
  
  candidates <- c("purity", "Purity", "tumorPurity", "tumourPurity")
  col <- candidates[candidates %in% names(pur)][1]
  if (is.na(col)) stop("Purity column not found. Columns: ", paste(names(pur), collapse = ", "))
  
  as.numeric(pur[[col]][1])
}

get_cnv_gr <- function(c) {
  # c can be file path OR data.frame
  cnv <- if (is.character(c) && length(c) == 1) {
    read.delim(c, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    c
  }
  
  # Try to auto-detect common PURPLE column names
  chr_col   <- c("chromosome","chr","Chromosome","contig")[c("chromosome","chr","Chromosome","contig") %in% names(cnv)][1]
  start_col <- c("start","startPos","Start","posStart")[c("start","startPos","Start","posStart") %in% names(cnv)][1]
  end_col   <- c("end","endPos","End","posEnd")[c("end","endPos","End","posEnd") %in% names(cnv)][1]
  
  maj_col   <- c("majorAlleleCopyNumber","major_cn","majorCopyNumber","majorAlleleCN")[c("majorAlleleCopyNumber","major_cn","majorCopyNumber","majorAlleleCN") %in% names(cnv)][1]
  min_col   <- c("minorAlleleCopyNumber","minor_cn","minorCopyNumber","minorAlleleCN")[c("minorAlleleCopyNumber","minor_cn","minorCopyNumber","minorAlleleCN") %in% names(cnv)][1]
  
  if (any(is.na(c(chr_col,start_col,end_col,maj_col,min_col)))) {
    stop("CNV columns not found. Columns: ", paste(names(cnv), collapse = ", "))
  }
  
  GRanges(
    seqnames = cnv[[chr_col]],
    ranges   = IRanges(start = cnv[[start_col]], end = cnv[[end_col]]),
    major_cn = as.numeric(cnv[[maj_col]]),
    minor_cn = as.numeric(cnv[[min_col]])
  )
}

assign_cn_to_ssm <- function(ssm, seg_gr) {
  snv_gr <- GRanges(
    seqnames = ssm$chr,
    ranges = IRanges(start = ssm$pos, end = ssm$pos)
  )
  
  hits <- findOverlaps(snv_gr, seg_gr)
  
  ssm$major_cn <- NA_real_
  ssm$minor_cn <- NA_real_
  
  ssm$major_cn[queryHits(hits)] <- mcols(seg_gr)$major_cn[subjectHits(hits)]
  ssm$minor_cn[queryHits(hits)] <- mcols(seg_gr)$minor_cn[subjectHits(hits)]
  
  subset(ssm, !is.na(major_cn) & !is.na(minor_cn))
}

run_ccube_all <- function(ssm_list, cnvs, purities, init_k = 3, clonal_thr = 0.9) {
  
  stopifnot(setequal(names(ssm_list), names(cnvs)))
  stopifnot(setequal(names(ssm_list), names(purities)))
  
  sids <- intersect(intersect(names(ssm_list), names(cnvs)), names(purities))
  
  results <- vector("list", length(sids))
  names(results) <- sids
  
  for (sid in sids) {
    
    ssm <- ssm_list[[sid]]
    purity_value <- get_purity_value(purities[[sid]])
    seg_gr <- get_cnv_gr(cnvs[[sid]])
    
    ssm$purity <- purity_value
    ssm$normal_cn <- 2
    
    # assign CN to SNVs
    ssm <- assign_cn_to_ssm(ssm, seg_gr)
    
    if (nrow(ssm) < 50) {
      warning("Very few SNVs after CN overlap for ", sid, ": ", nrow(ssm))
    }
    
    ssm_core <- within(ssm, {
      major_cn_sub1 <- major_cn
      minor_cn_sub1 <- minor_cn
      major_cn_sub2 <- -100
      minor_cn_sub2 <- -100
      frac_cn_sub1  <- 1
      frac_cn_sub2  <- 0
      subclonal_cn  <- FALSE
      total_cn      <- major_cn + minor_cn
    })
    
    res <- CcubeCore(ssm_core, init = init_k, fit_mult = TRUE, fit_hyper = TRUE, use="use_base", verbose = FALSE)
    
    ssm2 <- AnnotateCcubeResults(ssm_core, res)
    ssm2$is_clonal <- ssm2$ccube_ccf >= clonal_thr
    
    # per-sample summary (handy for reviewer response)
    summ <- data.frame(
      sample = sid,
      purity = purity_value,
      n_snvs = nrow(ssm2),
      frac_clonal = mean(ssm2$is_clonal),
      median_ccf = median(ssm2$ccube_ccf, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    results[[sid]] <- list(ssm = ssm2, res = res, summary = summ)
  }
  
  results
}

ccube_out <- run_ccube_all(ssm_list, cnvs, purities, init_k = 3, clonal_thr = 0.9)


make_mydata_m6 <- function(ssm, purity_value,
                           normal_cn = 2,
                           min_depth = 10,
                           round_cn = TRUE) {
  
  stopifnot(all(c("ref_counts","var_counts") %in% names(ssm)))
  
  # If major/minor CN not present, stop (you must overlap CN segments first)
  if (!all(c("major_cn","minor_cn") %in% names(ssm))) {
    stop("ssm must contain major_cn and minor_cn. Overlap SNVs with PURPLE CNV segments first.")
  }
  
  mydata <- data.frame(
    mutation_id = if ("mutation_id" %in% names(ssm)) ssm$mutation_id else NA_character_,
    chr = if ("chr" %in% names(ssm)) ssm$chr else NA_character_,
    pos = if ("pos" %in% names(ssm)) ssm$pos else NA_integer_,
    
    ref_counts = as.integer(ssm$ref_counts),
    var_counts = as.integer(ssm$var_counts),
    
    normal_cn  = as.numeric(normal_cn),
    major_cn   = as.numeric(ssm$major_cn),
    minor_cn   = as.numeric(ssm$minor_cn),
    purity     = as.numeric(purity_value),
    
    stringsAsFactors = FALSE
  )
  
  # QC filters that prevent NaNs in the likelihood
  mydata$total_counts <- mydata$ref_counts + mydata$var_counts
  mydata <- subset(mydata, !is.na(purity) & purity > 0 & purity < 1)
  mydata <- subset(mydata, !is.na(major_cn) & !is.na(minor_cn))
  mydata <- subset(mydata, major_cn >= 0 & minor_cn >= 0)
  mydata <- subset(mydata, total_counts >= min_depth)
  
  if (round_cn) {
    mydata$major_cn <- pmax(0, round(mydata$major_cn))
    mydata$minor_cn <- pmax(0, round(mydata$minor_cn))
  }
  
  # Remove loci with total tumour CN < 1 (these break the model)
  mydata$total_cn <- mydata$major_cn + mydata$minor_cn
  mydata <- subset(mydata, total_cn >= 1)
  
  # Drop helper columns if you want
  mydata$total_counts <- NULL
  mydata$total_cn <- NULL
  
  mydata
}

sid <- names(ssm_list)[1]
ssm <- ssm_list[[sid]]

purity_value <- get_purity_value(purities[[sid]])  # you already have this helper
seg_gr <- get_cnv_gr(cnvs[[sid]])
ssm <- assign_cn_to_ssm(ssm, seg_gr)               # adds major_cn/minor_cn

mydata <- make_mydata_m6(ssm, purity_value)

mydata_core <- within(mydata, {
  major_cn_sub1 <- major_cn
  minor_cn_sub1 <- minor_cn
  major_cn_sub2 <- -100
  minor_cn_sub2 <- -100
  frac_cn_sub1  <- 1
  frac_cn_sub2  <- 0
  subclonal_cn  <- FALSE
  total_cn      <- major_cn + minor_cn
})

res <- CcubeCore(mydata_core,
                 init = 3,
                 fit_mult = TRUE,
                 fit_hyper = TRUE,
                 use = "use_one",
                 verbose = TRUE)

ssm2 <- AnnotateCcubeResults(mydata_core, res)
ssm2$is_clonal <- ssm2$ccube_ccf >= 0.99


ccube_results <- list()
ssm_results <- list()

sample_ids <- names(ssm_list)

for (sid in sample_ids) {
  
  message("Running: ", sid)
  
  try({
    
    ## --- get data ---
    ssm <- ssm_list[[sid]]
    purity_value <- get_purity_value(purities[[sid]])
    seg_gr <- get_cnv_gr(cnvs[[sid]])
    
    ## --- assign CN ---
    ssm <- assign_cn_to_ssm(ssm, seg_gr)
    
    ## --- build input ---
    mydata <- make_mydata_m6(ssm, purity_value)
    
    ## --- convert to CcubeCore format ---
    mydata_core <- within(mydata, {
      major_cn_sub1 <- major_cn
      minor_cn_sub1 <- minor_cn
      major_cn_sub2 <- -100
      minor_cn_sub2 <- -100
      frac_cn_sub1  <- 1
      frac_cn_sub2  <- 0
      subclonal_cn  <- FALSE
      total_cn      <- major_cn + minor_cn
    })
    
    ## --- run ccube ---
    res <- CcubeCore(
      mydata_core,
      init = 3,
      fit_mult = TRUE,
      fit_hyper = TRUE,
      use = "use_one",
      verbose = FALSE
    )
    
    ## --- annotate results ---
    ssm2 <- AnnotateCcubeResults(mydata_core, res)
    ssm2$is_clonal <- ssm2$ccube_ccf = 1
    ssm2$sample_id <- sid
    ## --- store ---
    ccube_results[[sid]] <- res
    ssm_results[[sid]] <- ssm2
    
  }, silent = FALSE)
}

all_ssm <- do.call(rbind, ssm_results)
tapply(all_ssm$is_clonal, all_ssm$sample_id, sum, na.rm = TRUE)


clonal_counts <- tapply(all_ssm$is_clonal, all_ssm$sample_id, sum, na.rm = TRUE)

clonal_df <- data.frame(
  Sample = names(clonal_counts),
  n_clonal = as.vector(clonal_counts)
)

#write.csv(clonal_df, "clonal_counts_per_sample_ccube.csv", row.names = FALSE)


