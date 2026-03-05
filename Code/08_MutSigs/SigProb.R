################################################################################
# Manuscript: Single-cell whole-genome sequencing reveals convergent evolution in Burkitt lymphoma
# Description: Script to look at mutational signature probabilities per SBS driver
# Author: Alexander Steemers
################################################################################

# load libraries

library(dplyr)
library(GenomicRanges)
library(GenomeInfoDb)
library(MutationalPatterns)
library(BSgenome.Hsapiens.NCBI.GRCh38)

ref_genome <- BSgenome.Hsapiens.NCBI.GRCh38

drivers <- read.csv("mutation_data_filtered.csv",
                    stringsAsFactors = FALSE)  # generated after running Drivers_bulk.R script

## 1) Filter SNVs only
snv_data <- drivers %>%
  filter(!is.na(CHROM), !is.na(POS), !is.na(REF), !is.na(ALT)) %>%
  filter(nchar(REF) == 1, nchar(ALT) == 1) %>%
  filter(!Sample %in% c("P856_R2", "P856_R1", "PRN4_D2")) # remove samples for which other samples from the same patient already exist (they have the same set of drivers)

snv_data <- snv_data[!duplicated(snv_data[, c("CHROM", "POS", "REF", "ALT")]), ]

# clean types
snv_data$CHROM  <- as.character(snv_data$CHROM)
snv_data$POS    <- as.integer(snv_data$POS)
snv_data$REF    <- as.character(snv_data$REF)
snv_data$ALT    <- as.character(snv_data$ALT)
snv_data$Sample <- as.character(snv_data$Sample)

## 2) Build GRanges
gr <- GRanges(
  seqnames = snv_data$CHROM,
  ranges   = IRanges(start = snv_data$POS, width = 1),
  strand   = "*",
  REF      = snv_data$REF,
  ALT      = snv_data$ALT,
  Sample   = snv_data$Sample
)

# add gene if available (SYMBOL in your file)
if ("SYMBOL" %in% colnames(snv_data)) {
  mcols(gr)$gene <- as.character(snv_data$SYMBOL)
} else {
  mcols(gr)$gene <- NA_character_
}

## 3) Make each mutation a separate GRangesList element
grl_mut <- split(gr, seq_along(gr))  # each element has exactly 1 mutation

## Name each mutation (these names become mut_matrix column names)
names(grl_mut) <- paste0(
  ifelse(is.na(mcols(gr)$gene) | mcols(gr)$gene == "", "NA", mcols(gr)$gene), "_",
  mcols(gr)$Sample, "_",
  as.character(seqnames(gr)), ":",
  start(gr), "_",
  mcols(gr)$REF, ">", mcols(gr)$ALT
)

## 4) Force NCBI-style seqnames: remove "chr" prefix if present
grl_mut <- endoapply(grl_mut, function(x) {
  old <- seqlevels(x)
  new <- sub("^chr", "", old)
  names(new) <- old
  renameSeqlevels(x, new)
})

## 5) Set genome + keep only seqlevels present in reference
ref_seqlevels <- seqlevels(ref_genome)
grl_mut <- endoapply(grl_mut, function(x) {
  genome(x) <- genome(ref_genome)[1]  # "GRCh38"
  keepSeqlevels(x, intersect(seqlevels(x), ref_seqlevels), pruning.mode = "coarse")
})

## Drop any empty elements (just in case)
grl_mut <- grl_mut[elementNROWS(grl_mut) > 0]
stopifnot(length(grl_mut) > 0)

## 6) Build SBS96 matrix: one column per mutation (ALL samples combined)
mut_mat_per_mut <- as.matrix(mut_matrix(vcf_list = grl_mut, ref_genome = ref_genome))

## 7) ALSO build a separate SBS96-per-mutation matrix for each sample
mut_sample <- vapply(grl_mut, function(x) mcols(x)$Sample[1], character(1))
grl_mut_by_sample <- split(grl_mut, mut_sample)

mut_mat_per_mut_by_sample <- lapply(grl_mut_by_sample, function(grl_s) {
  as.matrix(mut_matrix(vcf_list = grl_s, ref_genome = ref_genome))
})

mut_mat_per_mut_by_sample <- mut_mat_per_mut_by_sample[
  !names(mut_mat_per_mut_by_sample) %in% c("PC1A", "P856_D")
] 

#  PC1A was not included in the mutational pattern analysis because of inflated mutation burden

mut_mat_per_mut_by_sample <- lapply(mut_mat_per_mut_by_sample, function(mat) {
  cn <- colnames(mat)
  for (old in names(pmcid_map)) {
    cn <- gsub(old, pmcid_map[[old]], cn, fixed = TRUE)
  }
  colnames(mat) <- cn
  mat
})

all_signatures = get_known_signatures()
sbsblood <- read.table("Data/sigfit_cosmic3_bloodsig_Aug2020.txt", sep = "\t", header = T)
sbsblood = as.matrix(sbsblood)
SBSblood <- as.numeric(sbsblood[,"Signature.Blood"])
signatures <- cbind(SBSblood, all_signatures)

contri_mat <- read.csv("bulk_signature_contributions.csv", row.names = 1) # generated after running MutationalPatterns_PMC_Machado_pcawg.R script
contri_mat <- contri_mat[
  , !colnames(contri_mat) %in% c("PRN4_D2","P856_R1", "P856_R2")
] # remove samples for which other samples from the same patient already exist (they have the same set of drivers)

# signatures you want
sigs <- c("SBS1","SBS7a", "SBS9", "SBS17b", "SBS18", "SBSblood")
sub.signatures <- signatures[, sigs, drop = FALSE]
rownames(sub.signatures) = rownames(mut_mat_per_mut_by_sample$P3G6)
# sanity checks: make sure sample names match between the objects
common_samples <- intersect(names(mut_mat_per_mut_by_sample), colnames(contri_mat))
if (length(common_samples) == 0) {
  stop("No overlapping sample names between mut_mat_per_mut_by_sample and contri_mat.")
}

# optional: restrict to shared samples only
mut_mat_per_mut_by_sample2 <- mut_mat_per_mut_by_sample[common_samples]
contri_mat2 <- contri_mat[, common_samples, drop = FALSE]

# result table
list.mut <- data.frame(
  sample = character(),
  mutation = character(),
  signature = character(),
  probability = numeric(),
  stringsAsFactors = FALSE
)

# Loop over samples
for (samp in names(mut_mat_per_mut_by_sample2)) {
  
  mm <- mut_mat_per_mut_by_sample2[[samp]]  # SBS96 x mutations for this sample
  
  # Loop over mutations (columns)
  for (mut in colnames(mm)) {
    
    # find the SBS96 context/bin for this mutation (row name where value == 1)
    bin <- rownames(mm)[which(mm[, mut] == 1)]
    
    # If something odd happens (e.g., no bin or multiple bins), skip
    if (length(bin) != 1) next
    
    # Loop over signatures
    for (sig in sigs) {
      
      # contribution of this signature in this sample
      contrib <- contri_mat2[sig, samp]
      
      # probability of this SBS96 bin under this signature
      sig_prob <- sub.signatures[bin, sig]
      
      # your per-mutation "prop"
      prop <- contrib * sig_prob
      
      list.mut <- rbind(
        list.mut,
        data.frame(
          sample = samp,
          mutation = mut,
          signature = sig,
          probability = prop,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

list.mut<-list.mut[!is.na(list.mut$probability), ]
list.mut$probability.norm <- prop.table(list.mut$probability)

prob_mat <- list.mut %>%
  dplyr::select(mutation, signature, probability) %>%
  group_by(mutation, signature) %>%
  summarise(probability = sum(probability), .groups = "drop") %>%   # in case duplicates
  pivot_wider(
    names_from  = signature,
    values_from = probability,
    values_fill = 0
  ) %>%
  column_to_rownames("mutation") %>%
  as.matrix()

rs <- rowSums(prob_mat)
prob_mat_norm <- prob_mat
prob_mat_norm[rs > 0, ] <- prob_mat[rs > 0, ] / rs[rs > 0]
rowSums(prob_mat_norm)

pdf("sig_contri_per_driver_heatmap.pdf", width = 8, height = 10)

pheatmap(
  prob_mat_norm,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "orange", "red"))(100),
  fontsize_col = 10,
  fontsize_row  = 5,
  breaks = seq(0, 1, length.out = 101),        # forces color scale 0→1
  legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  legend_labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")
)

dev.off()

# get the number of drivers per signature that have high probability (>70%) of being caused by that specific signature

colSums(prob_mat_norm > 0.70)

