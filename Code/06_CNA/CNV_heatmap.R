################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Script to generate CNV-heatmaps for all patients
# Author: Alexander Steemers
################################################################################

# Load required packages 
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library(wesanderson)
library(rtracklayer)
library(GWASTools) #BiocManager::install("GWASTools")
library(biomaRt)

# Set working directory and define export directories

workDir <- # define
inputDir <- # define
plotDir <- # define
fileDir <- # define

setwd(workDir)

# Load hg38 centromere coordinates 

# NB PTATO centromere coordinates located at:
# Data/centromeres_hg38.csv
readr::write_csv(centromeres.hg38)

hg38Centromeres <- read.csv(paste0(metaDir, "centromeres_hg38.csv"),
                            sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Load gene coordinates file 

# NB Gencode table via SCGF
# Data/hg38_gencode_v27.txt
genePositions <- read.table(file = hg38_gencode_v27.txt)
colnames(genePositions) <- c("Gene", "Chromosome", "Start", "End")
genePositions$Chromosome <- gsub("chr", "", genePositions$Chromosome)

# Switch to BioMart in order to filter on protein coding genes

# Generate de novo genomic positions file via biomaRt
# Ensembl release 98 corresponds to GENCODE v32
#listEnsemblArchives() #to view available versions
# Select BioMart database and dataset
#listDatasets(ensembl)
#searchDatasets(mart = ensembl, pattern = "hsapiens")
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 115)
# Define attributes
attributes <- c("external_gene_name", "chromosome_name", "start_position", "end_position", "gene_biotype")

# Perform the BioMart query
genePositions <- getBM(attributes = attributes,
                       #filters = "external_gene_name",
                       #values = geneList,
                       mart = ensembl)
genePositions <- genePositions[which(genePositions$gene_biotype == "protein_coding"), ]

colnames(genePositions) <- c("GeneName", "Chromosome", "Start", "End", "Biotype")

# 0) Shared configuration  

# Patients are configured at the bottom; everything else is shared.

chromosomeNames <- as.character(1:22)
chromosomeNum   <- length(chromosomeNames)

# Pericentromere exclusion (same for all chrs)
periCentromeres <- rep(10000000, chromosomeNum)

# Chromosome ends (your values)
chrEnds <- c(248956422, 242193529, 198295559, 190214555, 181500000,
             170800000, 159345973, 145100000, 138394717, 133797422,
             135086622, 133275309, 114364328, 106900000, 101991189,
             90300000, 83257441, 80300000, 58617616, 64400000,
             46700000, 50818468)

chrMidPoints <- chrEnds / 2
chrCum <- c(0, cumsum(chrEnds))
chrTot <- chrCum[length(chrCum)]
chrMidCum <- cumsum(chrEnds) - chrMidPoints

# Plot/QC settings
smoothBin <- 50  # must be even
sdBin <- smoothBin

interThreshold <- 0.275
failThreshold  <- 0.4

# CNV colours/classes
cnvClasses <- c("Loss", "Gain", "LOH")
cnvColours <- c("Steelblue", "#E84B3D", "lightgoldenrod1")

bafColour     <- "orangered1"
sdColour      <- "darkblue"
failThresCol  <- "darkgrey"

# File patterns 
# these are generated after running PTATO
filePattern <- "*.readcounts.filtered.100kb.txt"
namePattern <- ".readcounts.filtered.100kb.txt"
binSize     <- "100kb"

# 1) Centromere coordinates (once, global)  
# Assumes hg38Centromeres already exists in your environment
hg38Centromeres <- hg38Centromeres[hg38Centromeres$chrom %in% chromosomeNames, ]

centromereCoords <- hg38Centromeres
centromereCoords$left.base  <- hg38Centromeres$left.base  + chrCum[-(chromosomeNum + 1)]
centromereCoords$right.base <- hg38Centromeres$right.base + chrCum[-(chromosomeNum + 1)]

# Exclude acrocentric coords
acroIdx <- c(13, 14, 15, 21, 22)
centromereCoords$left.base[acroIdx] <- NA

# 2) Helpers: smoothing + per-chr work 

# Your readcount smoothing: centered mean over smoothBin/2 each side
smooth_readcount_centered <- function(x, smoothBin) {
  out <- rep(NA_real_, length(x))
  half <- smoothBin / 2
  if (length(x) < (smoothBin + 1)) return(out)
  for (i in (1 + half):(length(x) - half)) {
    out[i] <- mean(x[(i - half):(i + half)], na.rm = TRUE)
  }
  out
}

# Your BAF smoothing: keep your original behaviour (forward window)
smooth_baf_forward <- function(x, smoothBin) {
  out <- rep(NA_real_, length(x))
  half <- smoothBin / 2
  if (length(x) < (smoothBin + 1)) return(out)
  for (i in (1 + half):(length(x) - half)) {
    out[i] <- mean(x[i:(i + smoothBin)], na.rm = TRUE)
  }
  out
}

exclude_centromere_region <- function(df, chrCent, peri, chrCount, valueCol) {
  # Mask pericentromere
  inPeri <- (df$end.pos  > (chrCent$left.base  - peri[chrCount])) &
            (df$start.pos < (chrCent$right.base + peri[chrCount]))
  df[[valueCol]][inPeri] <- NA

  # Mask acrocentric short arms (as you had it)
  for (acro in c(13, 14, 15, 21, 22)) {
    df[[valueCol]][(df$Chromosome == acro) &
                   (df$end.pos < (chrCent$left.base - peri[chrCount]))] <- NA
  }
  df
}

add_cumulative_coords <- function(df, chrCount) {
  df$start.pos <- df$start.pos + chrCum[chrCount]
  df$end.pos   <- df$end.pos   + chrCum[chrCount]
  df
}

add_cumulative_coords_cnv <- function(df, chrCount) {
  df$start <- df$start + chrCum[chrCount]
  df$end   <- df$end   + chrCum[chrCount]
  df
}

run_patient_qc_plots <- function(patient_id,
                                 patient_dir,
                                 batch_dirs,
                                 bulk_name,
                                 metadata_csv,
                                 inputDir,
                                 metaDir,
                                 covYaxis,
                                 bafYaxis,
                                 axisStep,
                                 altBackground = FALSE) {
  
  safe_read_table <- function(path, header = TRUE) {
    if (!file.exists(path)) return(NULL)
    tryCatch(read.table(file = path, header = header), error = function(e) NULL)
  }
  
  safe_read_csv <- function(path) {
    if (!file.exists(path)) return(NULL)
    tryCatch(readr::read_csv(file = path, show_col_types = FALSE), error = function(e) NULL)
  }
  
  # Source once per patient run
  source("Data/generalFunctionsARB.R")
  source("Data/readFunctions.R")
  
  plotDir_patient <- paste0("ChromosomePlots/", patient_dir, "Plots/")
  dir.create(plotDir_patient, recursive = TRUE, showWarnings = FALSE)
  dir.create(fileDir_patient, recursive = TRUE, showWarnings = FALSE)
  
  # ---- 1) Read metadata (order source of truth) ----
  treeDF <- safe_read_csv(metadata_csv)
  
  # Try to autodetect sample column
  candidate_cols <- c("sample", "Sample", "sampleName", "SampleName", "cell", "Cell", "barcode", "Barcode")
  sample_col <- intersect(candidate_cols, colnames(treeDF))
  if (length(sample_col) == 0) {
    # fallback: pick first column that is character-like and has many unique values
    char_cols <- colnames(treeDF)[sapply(treeDF, function(x) is.character(x) || is.factor(x))]
    if (length(char_cols) == 0) stop("No character/factor columns in metadata CSV to use as sample list.")
    # choose the one with max unique
    uniq_counts <- sapply(treeDF[char_cols], function(x) length(unique(as.character(x))))
    sample_col <- char_cols[which.max(uniq_counts)]
  } else {
    sample_col <- sample_col[1]
  }
  
  ordered_samples <- as.character(treeDF[[sample_col]])
  ordered_samples <- ordered_samples[!is.na(ordered_samples) & ordered_samples != ""]
  
  if (length(ordered_samples) == 0) {
    stop("No sample names found in metadata column: ", sample_col)
  }
  
  message("Patient ", patient_id, ": using metadata column '", sample_col,
          "' with ", length(ordered_samples), " samples (row order preserved).")
  
  # ---- 2) Index all readcount files across batches once ----
  all_readcount_files <- unlist(lapply(batch_dirs, function(batch) {
    list.files(
      path = paste0(inputDir, patient_dir, batch),
      pattern = filePattern,
      full.names = TRUE,
      recursive = TRUE
    )
  }))
  all_readcount_files <- all_readcount_files[!grepl("old", all_readcount_files)]
  
  # Map file -> sampleName (same logic you used)
  file_sample_names <- gsub(paste0(".+/(.+)", namePattern), "\\1", all_readcount_files)
  readcount_map <- setNames(all_readcount_files, file_sample_names)
  
  # ---- 3) Create ONE PDF per patient ----
  pdf_out <- paste0(plotDir_patient, patient_id, "_chromosome_plot_QC_ALLSAMPLES.pdf")
  
  # We’ll size by number of ordered samples (even if some missing -> still reserve page per plotted sample)
  pdf(pdf_out, height = (5 * length(ordered_samples)), width = chromosomeNum * 1.5)
  par(mfrow = c(length(ordered_samples), 1))
  par(mar = c(5.5, 5.5, 5.5, 5.5))
  
  # ---- 4) Loop in metadata order ----
  for (sampleName in ordered_samples) {
    # Locate readcount file
    cell <- readcount_map[[sampleName]]
    
    if (is.null(cell) || is.na(cell) || !file.exists(cell)) {
      # placeholder page so ordering stays intact
      plot(0, 0, type = "n", xaxt = "n", yaxt = "n",
           xlab = "", ylab = "",
           main = paste0(sampleName, " (readcount file not found)"))
      next
    }
    
    message("  Plotting sample: ", sampleName)
    
    tryCatch({
      # Required: readcount table
      readCountTable <- safe_read_table(cell, header = TRUE)
      if (is.null(readCountTable)) {
        plot(0, 0, type = "n", xaxt = "n", yaxt = "n",
             xlab = "", ylab = "",
             main = paste0(sampleName, " (readcount unreadable)"))
        next
      }
      
      # Infer batch from file path (so CNV path matches the right batch)
      # Assumes ".../<patient_dir>/<batch>/...<sample>.readcounts..."
      batch_hit <- NA_character_
      for (b in batch_dirs) {
        if (grepl(paste0("/", b, "/"), cell, fixed = FALSE)) { batch_hit <- b; break }
      }
      if (is.na(batch_hit)) batch_hit <- batch_dirs[1]
      
      # Optional: CNV table (batch-specific)
      cnvFile <- paste0(inputDir, patient_dir, batch_hit,
                        "/intermediate/svs/Integration/", batch_hit, "/",
                        bulk_name, "/", sampleName, ".integrated.cnvs.txt")
      cnvTable <- safe_read_table(cnvFile, header = TRUE)
      if (!is.null(cnvTable)) {
        cnvTable <- cnvTable[cnvTable$chrom %in% chromosomeNames, ]
      }
      
      # Build chromosome-combined tables
      readCountTemp <- NULL
      cnvTableTemp  <- NULL
      
      for (chrCount in seq_along(chromosomeNames)) {
        chr <- chromosomeNames[chrCount]
        chrCent <- hg38Centromeres[hg38Centromeres$chrom == chr, ]
        
        rcChr <- readCountTable[readCountTable$Chromosome == chr, ]
        rcChr$medianReadCountSmooth <- smooth_readcount_centered(rcChr$medianReadCount, smoothBin)
        rcChr <- exclude_centromere_region(rcChr, chrCent, periCentromeres, chrCount, "medianReadCountSmooth")
        rcChr <- add_cumulative_coords(rcChr, chrCount)
        readCountTemp <- rbind(readCountTemp, rcChr)
        
        if (!is.null(cnvTable)) {
          cnvChr <- cnvTable[cnvTable$chrom == chr, ]
          if (nrow(cnvChr) > 0) {
            cnvChr <- add_cumulative_coords_cnv(cnvChr, chrCount)
            cnvTableTemp <- rbind(cnvTableTemp, cnvChr)
          }
        }
      }
      
      # Plot settings
      axisStep <- 100000000
      covYaxis <- c(0, 2)
      
      # Background line (white) to set up panel
      plot(readCountTemp$start.pos, readCountTemp$medianReadCountSmooth,
           lwd = 3, col = "white", type = "l",
           main = "", xlab = "", ylab = "",
           xaxt = "n", yaxt = "n",
           cex.axis = 2, cex.main = 4,
           xlim = c(0, chrTot))
      
      par(new = TRUE)
      
      # CNV regions
      if (!is.null(cnvTableTemp) && nrow(cnvTableTemp) > 0) {
        for (k in seq_along(cnvClasses)) {
          cnvSubset <- cnvTableTemp[cnvTableTemp$CopyNumber == cnvClasses[k], ]
          if (nrow(cnvSubset) > 0) {
            rect(xleft = cnvSubset$start, ybottom = -1,
                 xright = cnvSubset$end, ytop = 3,
                 col = cnvColours[k], border = NA)
          }
        }
      }
      
      par(new = TRUE)
      
      # FINAL overlay plot (no sample name, no axis labels/ticks)
      plot(readCountTemp$start.pos, readCountTemp$medianReadCountSmooth,
           lwd = 3, type = "l",
           main = "", xlab = "", ylab = "",
           xaxt = "n", yaxt = "n",
           cex.axis = 2, cex.main = 4,
           xlim = c(0, chrTot), ylim = covYaxis)
      
      # Keep vertical guides if you want
      abline(v = chrCum, lwd = 3)
      abline(v = centromereCoords$left.base  - periCentromeres, lty = 2, lwd = 3)
      abline(v = centromereCoords$right.base + periCentromeres, lty = 2, lwd = 3)
      
    }, error = function(e) {
      plot(0, 0, type = "n", xaxt = "n", yaxt = "n",
           xlab = "", ylab = "",
           main = paste0(sampleName, " (ERROR: ", conditionMessage(e), ")"))
    })
  }
  
  dev.off()
  message("Wrote: ", pdf_out)
}

# 4) Patient configs + run everything
# here the patients have different batches because the samples were split into batches when running PTATO
# adjust appropriately

patients <- list(
  list(id="P3G6", dir="P3G6/", batches=c("batch1","batch2"), bulk="PB11197-BLBM-MSCBULK",
       meta="Data/samples_treeposition_P3G6.csv"),
  list(id="PIA9", dir="PIA9/", batches=paste0("batch",1:7), bulk="PIA9GBDBMS79",
       meta="Data/samples_treeposition_PIA9.csv"),
  list(id="PJBU", dir="PJBU/", batches=paste0("batch",1:8), bulk="PJBUGBDBMS83",
       meta="Data/samples_treeposition_PJBU.csv"),
  list(id="PRN4", dir="PRN4/", batches=paste0("batch",1:3), bulk="PB08410-BM-MSCBULK",
       meta="Data/samples_treeposition_PRN4.csv"),
  list(id="P856", dir="P856/", batches=paste0("batch",1:4), bulk="PB14458-BLBM-MSCBULK",
       meta="Data/samples_treeposition_P856.csv"),
  list(id="PVA9", dir="PVA9/", batches=paste0("batch",1:8), bulk="PVA9GBDBMS79",
       meta="Data/samples_treeposition_PVA9.csv")
)

for (p in patients) {
  run_patient_qc_plots(
    patient_id   = p$id,
    patient_dir  = p$dir,
    batch_dirs   = p$batches,
    bulk_name    = p$bulk,
    metadata_csv = p$meta,
    inputDir     = inputDir,
    metaDir      = metaDir,
    covYaxis     = covYaxis,
    bafYaxis     = bafYaxis,
    axisStep     = axisStep
  )
}

