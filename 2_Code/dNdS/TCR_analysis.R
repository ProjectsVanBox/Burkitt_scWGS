################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to compare strand bias
# Author: Alexander Steemers
# Date: August 2025
################################################################################

# Load libraries

library(dplyr)
library(VariantAnnotation)
library(GenomeInfoDb)
library(BiocGenerics)
library(ggplot2)
library(dndscv)
library(purrr)
library(tidyr)
library(tibble)
library(readr)
library(ggrepel)
library(gridExtra)
library(MutationalPatterns)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# Load grl lists

branch_grl_P3G6 <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_grl_P3G6.rds")
branch_grl_PVA9 <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_grl_PVA9.rds")
branch_grl_PIA9 <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_grl_PIA9.rds")
branch_grl_PJBU <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_grl_PJBU.rds")
branch_grl_P856 <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_grl_P856.rds")
branch_grl_PRN4 <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_grl_PRN4.rds")


# Define branches per sample

# P3G6
group1_P3G6 <- c("Z")
group2_P3G6 <- c("T","R","Y","X","P")
group3_P3G6 <- c("J","I","H","N","M","Q","S","W","V","U")
group4_P3G6 <- c("B","A","D","F","G","a")
branch_to_group_P3G6 <- c(
  setNames(rep("BL-Trunk",        length(group1_P3G6)), group1_P3G6),
  setNames(rep("BL-Intermediate", length(group2_P3G6)), group2_P3G6),
  setNames(rep("BL-Private",      length(group3_P3G6)), group3_P3G6),
  setNames(rep("WT",              length(group4_P3G6)), group4_P3G6)
)

# PVA9
group1_PVA9 <- c("J3")
group2_PVA9 <- c("H3","G3","q2","N","v2","E3")
group3_PVA9 <- c("O2","N2","Q2","S2","H2","G2","J2","L2","V2","F2","C2","B2","A2","a2","Z2","z","e2","u","t","w","r","q","o","n","c","b","e","h","g","k","a","Z","l2","k2","n2","Q","P","O","T","W","V","A3","z2","x2","w2","D3","s2","r2","u2","M","L","I3")
group4_PVA9 <- c("L3","K3","N3","K","G","F","I","T3","S3","B","A")
branch_to_group_PVA9 <- c(
  setNames(rep("BL-Trunk",        length(group1_PVA9)), group1_PVA9),
  setNames(rep("BL-Intermediate", length(group2_PVA9)), group2_PVA9),
  setNames(rep("BL-Private",      length(group3_PVA9)), group3_PVA9),
  setNames(rep("WT",              length(group4_PVA9)), group4_PVA9)
)

# PIA9
group1_PIA9 <- c("a2")
group2_PIA9 <- c("d2","E3","C3","t2","s2","e2","G3","I3","J3")
group3_PIA9 <- c("L","K","N","P","I","H","G","T","W","V","Y","B","A","M3","E","D","b","H3","I2","k2","n2","j2","q2","h2","g2","f2","z2","y2","x2","v2","u2","D3","c2","b2")
group4_PIA9 <- c("s","r","u","u","w","p","o","z","n","m","E2","D2","I","I2","K2","M2","O2","k","R2","j","i","V2","X2","h","c","d")
branch_to_group_PIA9 <- c(
  setNames(rep("BL-Trunk",        length(group1_PIA9)), group1_PIA9),
  setNames(rep("BL-Intermediate", length(group2_PIA9)), group2_PIA9),
  setNames(rep("BL-Private",      length(group3_PIA9)), group3_PIA9),
  setNames(rep("WT",              length(group4_PIA9)), group4_PIA9)
)

# PJBU
group1_PJBU <- c("m2")
group2_PJBU <- c("T","l2","j2","T2","S","Q","L")
group3_PJBU <- c("k","j","m","i","p","h","s","f","e","d","w","c","Z","Y","W","V","G2","F2","E2","D2","K2","M2","O2","C2","B2","A2","U2","W2","b2","a2","Z2","f2","e2","Y2","k2","N","M","P","R","K","J","I")
group4_PJBU <- c("p2","o2","r2","t2","v2","n2","A3","z2","C3","G3","F3","I3","G","F","D","A")
branch_to_group_PJBU <- c(
  setNames(rep("BL-Trunk",        length(group1_PJBU)), group1_PJBU),
  setNames(rep("BL-Intermediate", length(group2_PJBU)), group2_PJBU),
  setNames(rep("BL-Private",      length(group3_PJBU)), group3_PJBU),
  setNames(rep("WT",              length(group4_PJBU)), group4_PJBU)
)

# P856
group1_P856 <- c("A2") # BL trunk
group2_P856 <- c("z", "q", "o") # BL intermediate BM
group3_P856 <- c("W", "V", "Y", "U", "b", "T", "e", "j", "i", "h", "g", "K", "J", "M", "O", "I", "H", "E", "D", "C", "p", "w", "v", "u", "s", "r") # BL private
group4_P856 <- c("A", "B") # WT
branch_to_group_P856 <- c(
  setNames(rep("BL-Trunk",        length(group1_P856)), group1_P856),
  setNames(rep("BL-Intermediate", length(group2_P856)), group2_P856),
  setNames(rep("BL-Private",      length(group3_P856)), group3_P856),
  setNames(rep("WT",              length(group4_P856)), group4_P856)
)

# PRN4
group1_PRN4 <- c("h") # BL trunk 
group2_PRN4 <- c("p", "q", "a", "o", "m") # BL intermediate 
group3_PRN4 <- c("G", "F", "E", "C", "B", "A", "v", "u", "x", "t", "M", "L", "O", "Q", "S", "V", "U", "Z", "Y", "j", "i", "l", "n") # BL private 
group4_PRN4 <- c("g", "b", "c") # WT
branch_to_group_PRN4 <- c(
  setNames(rep("BL-Trunk",        length(group1_PRN4)), group1_PRN4),
  setNames(rep("BL-Intermediate", length(group2_PRN4)), group2_PRN4),
  setNames(rep("BL-Private",      length(group3_PRN4)), group3_PRN4),
  setNames(rep("WT",              length(group4_PRN4)), group4_PRN4)
)

# Bundle your per-patient objects and exact branch→group maps
grl_by_sample <- list(
  P3G6 = branch_grl_P3G6,
  PVA9 = branch_grl_PVA9,
  PIA9 = branch_grl_PIA9,
  PJBU = branch_grl_PJBU,
  P856 = branch_grl_P856,
  PRN4 = branch_grl_PRN4
)

branch2group_by_sample <- list(
  P3G6 = branch_to_group_P3G6,
  PVA9 = branch_to_group_PVA9,
  PIA9 = branch_to_group_PIA9,
  PJBU = branch_to_group_PJBU,
  P856 = branch_to_group_P856,
  PRN4 = branch_to_group_PRN4
)

# Annotate each branch with sample/branch/group; drop empties and any unmapped branches
annotate_with_map <- function(grl, sample, b2g, drop_empty = TRUE) {
  if (drop_empty) grl <- grl[elementNROWS(grl) > 0]
  if (length(grl) == 0) return(GRangesList())
  
  # Keep only branches that exist in the mapping
  keep_names <- intersect(names(grl), names(b2g))
  dropped <- setdiff(names(grl), keep_names)
  if (length(dropped)) {
    warning("Dropped unmapped branches in ", sample, ": ", paste(dropped, collapse=", "))
  }
  grl <- grl[keep_names]
  
  imap(grl, function(gr, br) {
    n <- length(gr)
    mcols(gr)$sample <- rep(sample, n)
    mcols(gr)$branch <- rep(br, n)
    mcols(gr)$group  <- rep(unname(b2g[[br]]), n)
    gr
  }) |> GRangesList()
}

annotated_per_sample <- imap(
  grl_by_sample,
  ~annotate_with_map(.x, .y, branch2group_by_sample[[.y]], drop_empty = TRUE)
)

# Combine all branches from all patients
all_branches_grl <- do.call(c, unname(annotated_per_sample))  # each element = a branch

# Determine each branch's group
branch_group <- vapply(
  all_branches_grl,
  function(gr) unique(as.character(mcols(gr)$group))[1],
  character(1)
)

# Build one GRangesList with 4 elements: Trunk / Intermediate / Private / WT
group_levels <- c("BL-Trunk","BL-Intermediate","BL-Private","WT")

group_grs <- lapply(group_levels, function(g) {
  idx <- which(branch_group == g)
  if (length(idx) == 0) {
    # return empty GRanges if that group isn't present
    GRanges()
  } else {
    unlist(all_branches_grl[idx], use.names = FALSE)
  }
})
names(group_grs) <- c("Trunk","Intermediate","Private","WT")

grl_by_group <- GRangesList(group_grs)

genes_hg38 <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Convert all GRanges in grl_by_group to UCSC-style chromosome names (chr1, chr2, ...)
grl_by_group <- endoapply(grl_by_group, function(gr) {
  # Drop unused seqlevels first to avoid errors
  gr <- keepSeqlevels(gr, seqlevelsInUse(gr), pruning.mode = "coarse")
  # Change naming style to UCSC (adds "chr")
  seqlevelsStyle(gr) <- "UCSC"
  gr
})

mut_mat_s <- mut_matrix_stranded(grl_by_group, ref_genome, genes_hg38)
mut_mat_s[1:5, 1:4]

branch <- c("Trunk", "Intermediate","Private", "WT" )

strand_counts <- strand_occurrences(mut_mat_s, by = branch)

strand_bias <- strand_bias_test(strand_counts)

ps1 <- plot_strand(strand_counts, mode = "relative")

ps2 <- plot_strand_bias(strand_bias, sig_type = "p")

grid.arrange(ps1, ps2)


strand_bias_notstrict <- strand_bias_test(strand_counts,
                                          p_cutoffs = c(0.5, 0.1, 0.05),
                                          fdr_cutoffs = 0.5
)
plot_strand_bias(strand_bias_notstrict, sig_type = "p")


repli_file <- system.file("extdata/ReplicationDirectionRegions.bed",
                          package = "MutationalPatterns"
)
repli_strand <- read.table(repli_file, header = TRUE)
# Store in GRanges object
repli_strand_granges <- GRanges(
  seqnames = repli_strand$Chr,
  ranges = IRanges(
    start = repli_strand$Start + 1,
    end = repli_strand$Stop
  ),
  strand_info = repli_strand$Class
)
# UCSC seqlevelsstyle
seqlevelsStyle(repli_strand_granges) <- "UCSC"
repli_strand_granges

repli_strand_granges$strand_info <- factor(repli_strand_granges$strand_info,
                                           levels = c("right", "left")
)

mut_mat_s_rep <- mut_matrix_stranded(grl_by_group, ref_genome, repli_strand_granges,
                                     mode = "replication"
)
strand_counts_rep <- strand_occurrences(mut_mat_s_rep, by = branch)
strand_bias_rep <- strand_bias_test(strand_counts_rep)

ps1 <- plot_strand(strand_counts_rep, mode = "relative")
ps2 <- plot_strand_bias(strand_bias_rep)
grid.arrange(ps1, ps2)
