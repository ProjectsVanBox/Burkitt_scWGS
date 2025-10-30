################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script for plotting MiXCR results for patient PIA9
# Author: Alexander Steemers
# Date: June 2025
################################################################################

# Load libraries

library(dplyr)
library(ComplexHeatmap)
library(randomcoloR)
library(circlize)
library(readxl)

# Set working directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MIXCR/")

# Load metadata

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx')

# read in MiXCR ouput for each patient

PIA9_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/MIXCR_new/PIA9/Output/.", pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)

# Read each .tsv file into a list

tsv_files_list <- lapply(PIA9_files, function(x) read.delim(x, header = TRUE, sep = "\t"))

# Set names of the list based on file names (without .tsv extension)

names(tsv_files_list) <- gsub("\\.tsv$", "", basename(PIA9_files))

# Combine tsv file list into one dataframe

combined_df <- do.call("rbind", lapply(tsv_files_list, as.data.frame)) 
combined_df$SampleId <- rep(names(tsv_files_list), sapply(tsv_files_list, nrow))

# Extract sample id from tsv file names

combined_df$NovogeneName <- sapply(1:length(combined_df$SampleId), FUN=function(x) strsplit(combined_df$SampleId, split="_")[[x]][1])

# Extract chain from tsv file names

combined_df$Chain <- sapply(1:length(combined_df$SampleId), FUN=function(x) strsplit(combined_df$SampleId, split="_")[[x]][3])

# Match encrypted file names to legible cell names

combined_df$SampleName <- combined_df$NovogeneName

# Add myc translocation info

combined_df$Myc_translocation <- input_df$Myc_translocation_IGV[match(combined_df$NovogeneName, input_df$Sample_name)]

# Extract V, D and J genes

combined_df$V_gene <- sapply(1:nrow(combined_df), FUN=function(x) strsplit(combined_df$allVHitsWithScore, split="\\*00")[[x]][1])
combined_df$D_gene <- sapply(1:nrow(combined_df), FUN=function(x) strsplit(combined_df$allDHitsWithScore, split="\\*00")[[x]][1])
combined_df$J_gene <- sapply(1:nrow(combined_df), FUN=function(x) strsplit(combined_df$allJHitsWithScore, split="\\*00")[[x]][1])

# Merge genes to receptor recombination names

combined_df <- mutate(combined_df, BCR = paste(V_gene, D_gene, J_gene, sep = "_"))

# Add info on whether contig is functional (in frame and no stop codon)

combined_df <- combined_df %>% 
  mutate(functional = case_when(grepl("_", aaSeqCDR3) ~ "No",
                                grepl("\\*", aaSeqCDR3) ~ "No", 
                                TRUE ~ "Yes"))

# Save file

write.table(combined_df, "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MIXCR/Data/PIA9_combined_df.txt")

# Possibility to continue from saved df

combined_df <- read.table("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MIXCR/Data/PIA9_combined_df.txt")

# Plot 1: IgH

# Only functional heavy chains in single cells

IGH_chains <- combined_df %>% filter(Chain == "IGH")
IGH_chains <- IGH_chains[!grepl("No", IGH_chains$functional),]

# Heavy chain recombinations

heavy_output <- c(sort(unique(IGH_chains$V_gene)), sort(unique(IGH_chains$D_gene)), sort(unique(IGH_chains$J_gene)), sort(unique(IGH_chains$aaSeqCDR3)))

heavy_output <- heavy_output[!is.na(heavy_output)]

# Vector of all IGH genes and CDR3 AA sequences

heavy_output <- c(sort(unique(IGH_chains$BCR)), sort(unique(IGH_chains$aaSeqCDR3)))

# Loop through all cells

cell_names <- IGH_chains$SampleName

# Matrix of cell names

igh_mat <- matrix(data=NA, nrow=length(cell_names), ncol=length(heavy_output))

n_bcr_igh <- length(unique(IGH_chains$BCR))
n_cdr3_igh <- length(unique(IGH_chains$aaSeqCDR3))

# Loop through IGH genes

for(i in 1:length(cell_names)) {
  bcr <-  IGH_chains$BCR[IGH_chains$SampleName == cell_names[i]]
  if(length(bcr)==0) {
    igh_mat[i,c(1:n_bcr_igh)] <- 0
  } else {
    for(j in 1:length(unique(bcr))) {
      igh_mat[i, which(heavy_output == bcr[j])] <- 1
    }
  }
}

# Loop through CDR3 aa

for(i in 1:length(cell_names)) {
  cdr3 <-  IGH_chains$aaSeqCDR3[IGH_chains$SampleName == cell_names[i]]
  if(length(cdr3)==0) {
    igh_mat[i,c(n_bcr_igh+1:n_cdr3_igh)] <- 0
  } else {
    for(j in 1:length(unique(cdr3))) {
      igh_mat[i, which(heavy_output == cdr3[j])] <- 1
    }
  }
}

rownames(igh_mat) <- cell_names
colnames(igh_mat) <- heavy_output

# Make second matrix converting 0 and 1 to unique numbers for heatmap colors

igh_mat_cols <- igh_mat

for(i in 1:ncol(igh_mat_cols)) {
  for(j in 1:nrow(igh_mat_cols)) {
    if(!is.na(igh_mat_cols[j,i]) & igh_mat_cols[j,i] != 0) {
      igh_mat_cols[j,i] <- igh_mat_cols[j,i] + i 
    }
  }
}


# Pick maximally discrete colors for different BCRs

cols <- distinctColorPalette(ncol(igh_mat_cols), runTsne = FALSE)
col_discrete <- colorRamp2(c(0:(length(cols)+1)), c("grey90", "white", cols))

# Remove duplicate rows

igh_mat_cols <- igh_mat_cols[!duplicated(rownames(igh_mat_cols)), ]

# Label and order samples based on their Myc translocation

myc_translocation_neg <- combined_df %>%
  filter(Myc_translocation == "No") %>%
  distinct(SampleName) %>%  
  pull(SampleName)

myc_translocation_pos <- combined_df %>%
  filter(Myc_translocation == "Yes") %>%
  distinct(SampleName) %>%  
  pull(SampleName)

# Combine them in the desired order

ordered_samples <- c(myc_translocation_neg, myc_translocation_pos)

# Keep only samples which are in igh_mat_cols

ordered_samples_in_mat <- ordered_samples[ordered_samples %in% rownames(igh_mat_cols)]

# Reorder the data frame

igh_mat_cols <- igh_mat_cols[ordered_samples_in_mat, ]

# Give gorup names

sample_groups <- ifelse(rownames(igh_mat_cols) %in% c(myc_translocation_neg), "MYC::IGH-",
                        ifelse(rownames(igh_mat_cols) %in% myc_translocation_pos, "MYC::IGH+", "Missing data"))

# Define colors
group_colors <- c(
  "MYC::IGH-" = "#E69F00",
  "MYC::IGH+"= "#D55E00"
)

top_ha = HeatmapAnnotation(
  Sample = sample_groups,
  col = list(Sample = group_colors),
  annotation_name_side = "left"
)

# Plot heatmap for IgH

Heatmap(t(igh_mat_cols),
        na_col="white",
        col=col_discrete,
        color_space = "RBG",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        border=TRUE,
        column_title= "PIA9 IGH",
        row_names_side = "left", 
        top_annotation = top_ha,
        show_heatmap_legend = F,
        split=c(rep("IgH genes", n_bcr_igh), rep("CDR3aa", n_cdr3_igh)))


# Plot 2: IgK

# Only functional kappa light chains in single cells

IGK_chains <- combined_df %>% filter(Chain == "IGK")
IGK_chains <- IGK_chains[!grepl("No", IGK_chains$functional),]

# Kappa chain recombinations

kappa_output <- c(sort(unique(IGK_chains$V_gene)), sort(unique(IGK_chains$D_gene)), sort(unique(IGK_chains$J_gene)), sort(unique(IGK_chains$aaSeqCDR3)))

kappa_output <- kappa_output[!is.na(kappa_output)]

# Vector of all IGK genes and CDR3 AA sequences

kappa_output <- c(sort(unique(IGK_chains$BCR)), sort(unique(IGK_chains$aaSeqCDR3)))

# Loop through all cells

cell_names <- IGK_chains$SampleName

# Matrix of cell names

igk_mat <- matrix(data=NA, nrow=length(cell_names), ncol=length(kappa_output))

n_bcr_igk <- length(unique(IGK_chains$BCR))
n_cdr3_igk <- length(unique(IGK_chains$aaSeqCDR3))

# Loop through IGK genes

for(i in 1:length(cell_names)) {
  bcr <-  IGK_chains$BCR[IGK_chains$SampleName == cell_names[i]]
  if(length(bcr)==0) {
    igk_mat[i,c(1:n_bcr_igk)] <- 0
  } else {
    for(j in 1:length(unique(bcr))) {
      igk_mat[i, which(kappa_output == bcr[j])] <- 1
    }
  }
}

# Loop through CDR3 aa

for(i in 1:length(cell_names)) {
  cdr3 <-  IGK_chains$aaSeqCDR3[IGK_chains$SampleName == cell_names[i]]
  if(length(cdr3)==0) {
    igk_mat[i,c(n_bcr_igk+1:n_cdr3_igk)] <- 0
  } else {
    for(j in 1:length(unique(cdr3))) {
      igk_mat[i, which(kappa_output == cdr3[j])] <- 1
    }
  }
}

rownames(igk_mat) <- cell_names
colnames(igk_mat) <- kappa_output

# Make second matrix converting 0 and 1 to unique numbers for heatmap colors

igk_mat_cols <- igk_mat

for(i in 1:ncol(igk_mat_cols)) {
  for(j in 1:nrow(igk_mat_cols)) {
    if(!is.na(igk_mat_cols[j,i]) & igk_mat_cols[j,i] != 0) {
      igk_mat_cols[j,i] <- igk_mat_cols[j,i] + i 
    }
  }
}


# Pick maximally discrete colors for different BCRs

cols <- distinctColorPalette(ncol(igk_mat_cols), runTsne = FALSE)
col_discrete <- colorRamp2(c(0:(length(cols)+1)), c("grey90", "white", cols))

# Remove duplicate rows

igk_mat_cols <- igk_mat_cols[!duplicated(rownames(igk_mat_cols)), ]

# Keep only samples which are in igk_mat_cols

ordered_samples_in_mat <- ordered_samples[ordered_samples %in% rownames(igk_mat_cols)]

# Reorder the data frame

igk_mat_cols <- igk_mat_cols[ordered_samples_in_mat, ]

# Give gorup names

sample_groups <- ifelse(rownames(igk_mat_cols) %in% c(myc_translocation_neg), "MYC::IGH-",
                        ifelse(rownames(igk_mat_cols) %in% myc_translocation_pos, "MYC::IGH+", "Missing data"))

# Define colors
group_colors <- c(
  "MYC::IGH-" = "#E69F00",
  "MYC::IGH+"= "#D55E00"
)

top_ha = HeatmapAnnotation(
  Sample = sample_groups,
  col = list(Sample = group_colors),
  annotation_name_side = "left"
)

# Plot heatmap for IgK

Heatmap(t(igk_mat_cols),
        na_col="white",
        col=col_discrete,
        color_space = "RBG",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        border=TRUE,
        column_title= "PIA9 IGK",
        row_names_side = "left", 
        top_annotation = top_ha,
        show_heatmap_legend = F,
        split=c(rep("IgK genes", n_bcr_igk), rep("CDR3aa", n_cdr3_igk)))


# Plot 3: IgL

# Only functional lambda light chains in single cells

IGL_chains <- combined_df %>% filter(Chain == "IGL")
IGL_chains <- IGL_chains[!grepl("No", IGL_chains$functional),]

# Kappa chain recombinations

lambda_output <- c(sort(unique(IGL_chains$V_gene)), sort(unique(IGL_chains$D_gene)), sort(unique(IGL_chains$J_gene)), sort(unique(IGL_chains$aaSeqCDR3)))

lambda_output <- lambda_output[!is.na(lambda_output)]

# Vector of all IGL genes and CDR3 AA sequences

lambda_output <- c(sort(unique(IGL_chains$BCR)), sort(unique(IGL_chains$aaSeqCDR3)))

# Loop through all cells

cell_names <- IGL_chains$SampleName

# Matrix of cell names

igl_mat <- matrix(data=NA, nrow=length(cell_names), ncol=length(lambda_output))

n_bcr_igl <- length(unique(IGL_chains$BCR))
n_cdr3_igl<- length(unique(IGL_chains$aaSeqCDR3))

# Loop through IGL genes

for(i in 1:length(cell_names)) {
  bcr <-  IGL_chains$BCR[IGL_chains$SampleName == cell_names[i]]
  if(length(bcr)==0) {
    igl_mat[i,c(1:n_bcr_igl)] <- 0
  } else {
    for(j in 1:length(unique(bcr))) {
      igl_mat[i, which(lambda_output == bcr[j])] <- 1
    }
  }
}

# Loop through CDR3 aa

for(i in 1:length(cell_names)) {
  cdr3 <-  IGL_chains$aaSeqCDR3[IGL_chains$SampleName == cell_names[i]]
  if(length(cdr3)==0) {
    igl_mat[i,c(n_bcr_igl+1:n_cdr3_igl)] <- 0
  } else {
    for(j in 1:length(unique(cdr3))) {
      igl_mat[i, which(lambda_output == cdr3[j])] <- 1
    }
  }
}

rownames(igl_mat) <- cell_names
colnames(igl_mat) <- lambda_output

# Make second matrix converting 0 and 1 to unique numbers for heatmap colors

igl_mat_cols <- igl_mat

for(i in 1:ncol(igl_mat_cols)) {
  for(j in 1:nrow(igl_mat_cols)) {
    if(!is.na(igl_mat_cols[j,i]) & igl_mat_cols[j,i] != 0) {
      igl_mat_cols[j,i] <- igl_mat_cols[j,i] + i 
    }
  }
}


# Pick maximally discrete colors for different BCRs

cols <- distinctColorPalette(ncol(igl_mat_cols), runTsne = FALSE)
col_discrete <- colorRamp2(c(0:(length(cols)+1)), c("grey90", "white", cols))

# Remove duplicate rows

igl_mat_cols <- igl_mat_cols[!duplicated(rownames(igl_mat_cols)), ]

# Keep only samples which are in igl_mat_cols

ordered_samples_in_mat <- ordered_samples[ordered_samples %in% rownames(igl_mat_cols)]

# Reorder the data frame

igl_mat_cols <- igl_mat_cols[ordered_samples_in_mat, ]

# Give group names

sample_groups <- ifelse(rownames(igl_mat_cols) %in% c(myc_translocation_neg), "MYC::IGH-",
                        ifelse(rownames(igl_mat_cols) %in% myc_translocation_pos, "MYC::IGH+", "Missing data"))

# Define colors
group_colors <- c(
  "MYC::IGH-" = "#E69F00",
  "MYC::IGH+"= "#D55E00"
)

top_ha = HeatmapAnnotation(
  Sample = sample_groups,
  col = list(Sample = group_colors),
  annotation_name_side = "left"
)

# Plot heatmap for IGL

Heatmap(t(igl_mat_cols),
        na_col="white",
        col=col_discrete,
        color_space = "RBG",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        border=TRUE,
        column_title= "PIA9 IGL",
        row_names_side = "left", 
        top_annotation = top_ha,
        show_heatmap_legend = F,
        split=c(rep("IGL genes", 7), rep("CDR3aa", 7)))

# Plot heatmap for all Ig together

# Get the union of all unique rownames

all_cells <- unique(c(rownames(igh_mat_cols), rownames(igk_mat_cols), rownames(igl_mat_cols)))

# Create a full empty matrix with NAs 

igh_mat_aligned <- matrix(NA, nrow = length(all_cells), ncol = ncol(igh_mat_cols))
rownames(igh_mat_aligned) <- all_cells
colnames(igh_mat_aligned) <- colnames(igh_mat_cols)
matching_igh <- match(rownames(igh_mat_cols), all_cells)
igh_mat_aligned[matching_igh, ] <- igh_mat_cols

igk_mat_aligned <- matrix(NA, nrow = length(all_cells), ncol = ncol(igk_mat_cols))
rownames(igk_mat_aligned) <- all_cells
colnames(igk_mat_aligned) <- colnames(igk_mat_cols)
matching_igk <- match(rownames(igk_mat_cols), all_cells)
igk_mat_aligned[matching_igk, ] <- igk_mat_cols

igl_mat_aligned <- matrix(NA, nrow = length(all_cells), ncol = ncol(igl_mat_cols))
rownames(igl_mat_aligned) <- all_cells
colnames(igl_mat_aligned) <- colnames(igl_mat_cols)
matching_igl <- match(rownames(igl_mat_cols), all_cells)
igl_mat_aligned[matching_igl, ] <- igl_mat_cols

# Now combine them

bcr_mat <- cbind(igh_mat_aligned, igk_mat_aligned, igl_mat_aligned)

# Add colours

cols <- distinctColorPalette(ncol(igh_mat_cols)+ncol(igk_mat_cols)+ncol(igl_mat_cols), runTsne = FALSE)
col_discrete <- colorRamp2(c(0:(length(cols)+1)), c("grey90", "white", cols))

# Add missing cells 

ordered_samples <- as.character(ordered_samples)
rownames(bcr_mat) <- as.character(rownames(bcr_mat))
missing <- setdiff(ordered_samples, rownames(bcr_mat))

if (length(missing) > 0) {
  extra <- matrix(
    NA,                                 
    nrow = length(missing),
    ncol = ncol(bcr_mat),
    dimnames = list(missing, colnames(bcr_mat))
  )
  bcr_mat <- rbind(bcr_mat, extra)
}

# Order samples

bcr_mat <- bcr_mat[ordered_samples, ]


# Give group names

sample_groups <- ifelse(rownames(bcr_mat) %in% c(myc_translocation_neg), "MYC::IGH-",
                        ifelse(rownames(bcr_mat) %in% myc_translocation_pos, "MYC::IGH+", "Missing data"))

# Define colors

group_colors <- c(
  "MYC::IGH-" = "#E69F00",
  "MYC::IGH+"= "#D55E00"
)

top_ha = HeatmapAnnotation(
  Sample = sample_groups,
  col = list(Sample = group_colors),
  annotation_name_side = "left"
)

column_names_color <- group_colors[sample_groups]

# Plot heatmap for all BCRs

pdf("Figures/PIA9_MiXCR_BCR_complete_heatmap.pdf", width=8, height=14)
Heatmap(t(bcr_mat),
        na_col="white",
        col=col_discrete,
        color_space = "RBG",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        column_names_gp = gpar(col = column_names_color, fontsize = 3),
        border=TRUE,
        column_title= "PIA9 BCR recombinations",
        row_names_side = "right", 
        row_names_gp = gpar(fontsize = 6),
        show_heatmap_legend = F,
        row_title_rot = 0,
        top_annotation = top_ha,
        split=c(rep("IgH", n_bcr_igh), rep("IgH CDR3aa", n_cdr3_igh),
                rep("IgK", n_bcr_igk), rep("IgK CDR3aa", n_cdr3_igk),
                rep("IgL", n_bcr_igl), rep("IgL CDR3aa", n_cdr3_igl)))
dev.off()

bcr_mat_receptors_only <- bcr_mat[, grepl("^I", colnames(bcr_mat))]

pdf("Figures/PIA9_MiXCR_BCR_heatmap.pdf", width=8, height=14)
Heatmap(t(bcr_mat_receptors_only),
        na_col="white",
        col=col_discrete,
        color_space = "RBG",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = TRUE,
        column_names_gp = gpar(col = column_names_color, fontsize = 3),
        border=TRUE,
        column_title= "PIA9 BCR recombinations",
        row_names_side = "right", 
        show_heatmap_legend = F,
        row_title_rot = 0,
        top_annotation = top_ha,
        split=c(rep("IgH", n_bcr_igh),
                rep("IgK", n_bcr_igk),
                rep("IgL", n_bcr_igl)))
dev.off()

pdf("Figures/PIA9_MiXCR_BCR_heatmap.pdf_small.pdf", width = 3.5, height = 3.5)

# Reduce default margins
par(mar = c(2, 2, 2, 2))

Heatmap(
  t(bcr_mat_receptors_only),
  na_col = "white",
  col = col_discrete,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 4, col = column_names_color),  # smaller column labels
  row_names_gp = gpar(fontsize = 4),                                # smaller row labels
  border = TRUE,
  column_title = "PIA9 BCR recombinations",
  column_title_gp = gpar(fontsize = 7),                            # slightly reduced title font
  row_names_side = "right",
  show_heatmap_legend = FALSE,
  top_annotation = top_ha,
  heatmap_width = unit(3, "in"),
  heatmap_height = unit(3, "in")
)

dev.off()


