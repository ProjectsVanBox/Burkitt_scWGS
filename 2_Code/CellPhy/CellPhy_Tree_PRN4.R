################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at CellPhy Tree of Donor PRN4 (branch VCFs for manual checking, choosing percentage, signatures per branch)
# Author: Alexander Steemers
# Date: July 2025
################################################################################

# Load libraries and functions

library(tidyverse)
library(vcfR)
library(readxl)
library(ggtree)
library(treeio)
library(RColorBrewer)
library(stringi)
library(reshape2)
library(stringdist)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(ggnewscale)
ref_genome = 'BSgenome.Hsapiens.NCBI.GRCh38'
library(MutationalPatterns)
library(VariantAnnotation)
library(ggplot2)
library(dplyr)
library(cellPhyWrapperPlotting) #devtools::install_local("~/surfdrive/Shared/pmc_vanboxtel/general/2_Bioinformatics/Scripts/cellPhyWrapperPlotting/",force = TRUE)
source("~/surfdrive/Shared/pmc_vanboxtel/personal/asteemers/R_functions/plot_signature_contribution_alex.R")
source("~/surfdrive/Shared/pmc_vanboxtel/personal/asteemers/R_functions/plot_gg_tree_alex.R")
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')
theme_set(theme_classic())
source("~/surfdrive/Shared/pmc_vanboxtel/general/2_Bioinformatics/colors/Jurrians_colors.R")

# Load tree object

tree = readRDS("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PRN4/Filtered_samples/CPW_04/TreeObject0.4.RDS")
vcf = VariantAnnotation::readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PRN4/PRN4.vep.sub.SNV.autosomal.noBulksNoBlacklist.vcf")
vcf_bulk_LN = readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PRN4/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PRN4_bulk.vep.SMuRF.filtered.sorted.vcf.gz")
vcf_bulk_BM = readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/wgs/PMABM000FZO_PMABM000FZS_PMCRZ929TYS_WGS_PASS.vcf.gz")
vcf_bulk_BM_snv <- vcf_bulk_BM[ isSNV(vcf_bulk_BM)]
vcf_bulk_BM_snv_autosomal <- keepSeqlevels(vcf_bulk_BM_snv, paste0("chr", 1:22), pruning.mode = "coarse")

# prepare tree
tree = prepare_tree(tree)

# plot bare tree
p1 <- plot_gg_tree_base(tree)

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_absolute_tree_with_branch_info.pdf",
  plot = p1,
  width = 8,
  height = 6
)

p2 <- plot_gg_tree(tree, add_branch_length = F, add_bootstrap = F, add_tip_label = F,add_title = "") +
  geom_tree(color = "grey") +  # Set branch color to grey
  scale_x_continuous(name = "Mutation burden (SNVs)",
                     limits = c(0, 2600),
                     breaks = seq(0, 2600, by = 200)) +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 6)),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5, margin = margin(t = 6) )
  )

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_absolute_tree.pdf",
  plot = p2,
  width = 8,
  height = 6
)

p2b <- plot_gg_tree_alex(tree ,add_branch_length = F, add_bootstrap = T, add_tip_label = F,add_title = "", add_x = FALSE) +
  geom_tree(color = "grey") +  # Set branch color to grey
  scale_x_continuous(name = "Mutation burden (SNVs)",
                     limits = c(0, 2600),
                     breaks = seq(0, 2600, by = 200)) +
  theme_tree2() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 6)),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5, margin = margin(t = 6) )
  )

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_absolute_tree_with_bootstraps.pdf",
  plot = p2b,
  width = 8,
  height = 6
)

# Define functions to make callable loci-corrected trees

# Function to calculate the product of (1 - sensitivity) for each sample in a node
calculate_product_sensitivity <- function(samples, df_sensitivity) {
  # Split the 'samples' string by '|'
  sample_list <- strsplit(samples, "\\|")[[1]]
  
  # Retrieve sensitivities for these samples
  sensitivities <- df_sensitivity$sensitivity[df_sensitivity$Sample_name %in% sample_list]
  
  # Calculate the product of (1 - sensitivity)
  1 - prod(1 - sensitivities)
}

# function to take a tree and correct the branch lengths based on the product sensitivity
## by calculating the product of (1 - sensitivity) for each sample in a node
# requires samples (=samples column from output of cellphywrapper object: tree@data$samples .
# requires a CALLABLE df (=containing columns SAMPLES and CALLABLE (in nr of bases))

correct_branches <- function(tree, callable_df){
  # get sensitivity df by taking fraction of total possible callable loci
  max_callable <- 2745186691
  callable_df$sensitivity <- callable_df$Callable_Loci / max_callable
  sensitivity_df <- callable_df[c('Sample_name','sensitivity')]
  
  print(sensitivity_df)
  
  # calculate the sensitivity per node/branch using the calculate_product_sensitivity function
  tree@data$product_sensitivity <- sapply(tree@data$samples, calculate_product_sensitivity, sensitivity_df)
  
  # correct the branch lengths
  tree@data$corr_branch_lengths <- tree@data$branch_length / tree@data$product_sensitivity
  
  print(tree@data$branch_length)
  print(tree@data$corr_branch_lengths)
  
  tree@data$corr_branch_lengths[is.na(tree@data$corr_branch_lengths)] <- 0 # for merged germline sample
  
  # store also in the phylo object, required for downstream analyses based on the phylo object
  tree@phylo$edge.length <- tree@data[order(tree@data$node),]$corr_branch_lengths
  
  # store in branch_length column for cellphyplotting use
  tree@data$branch_length <- tree@data$corr_branch_lengths
  
  print(tree@phylo$edge.length)
  
  return(tree)
}

# Import callable loci and SBS PASS for donor

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') #dataframe
input_df_PRN4 <- input_df[input_df$Novogene_ID == "PRN4" & input_df$ResolveDNA_version %in% c("v1", "v2.0", "v2"),]
below_curve_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv")
low_call_frac_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv")
fail_vaf_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt")
filtered_samples <- unique(c(low_call_frac_df$Sample_name, below_curve_df$Sample_name, fail_vaf_df$samplename))  # samples that didn't pass QC
input_df_PRN4_filtered <- input_df_PRN4 %>% filter(!Sample_name %in% filtered_samples)
input_df_PRN4_filtered$Callable_Loci <- as.numeric(input_df_PRN4_filtered$Callable_Loci)

# Correct the branches
corrected_tree <- correct_branches(tree, input_df_PRN4_filtered)
corrected_tree@data$branch_length <- round(corrected_tree@data$branch_length)

# Plot corrected tree

p3 <- plot_gg_tree(corrected_tree, add_branch_length = F, add_bootstrap = F, add_tip_label = F,add_title = "") +
  geom_tree(color = "grey") +  # Set branch color to grey
  scale_x_continuous(name = "Mutation burden (SNVs)",
                     limits = c(0, 2600),
                     breaks = seq(0, 2600, by = 200)) +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 6)),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5, margin = margin(t = 6) )
  )
# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_absolute_tree_corrected_for_callable_loci.pdf",
  plot = p3,
  width = 8,
  height = 6
)

p3b <- plot_gg_tree_alex(corrected_tree ,add_branch_length = F, add_bootstrap = T, add_tip_label = F,add_title = "", add_x = FALSE) +
  geom_tree(color = "grey") +  # Set branch color to grey
  scale_x_continuous(name = "Mutation burden (SNVs)",
                     limits = c(0, 2600),
                     breaks = seq(0, 2600, by = 200)) +
  theme_tree2() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 6)),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5, margin = margin(t = 6) )
  )

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_absolute_tree_with_bootstraps_corrected_for_callable_loci.pdf",
  plot = p3b,
  width = 8,
  height = 6
)

x <- as_tibble(corrected_tree)

# 1. Get all non-NA tip labels
tip_labels <- x$tip.label[!is.na(x$tip.label)]

# 2. For each tip label, sum the branch lengths where it's mentioned in `samples`
mutation_df <- tibble(TipLabel = tip_labels) %>%
  rowwise() %>%
  mutate(
    TotalMutations = sum(
      x$corr_branch_lengths[
        str_detect(x$samples, stringr::fixed(TipLabel))
      ],
      na.rm = TRUE
    )
  ) %>%
  ungroup()


# VCF per branch
branch_vcf = extract_vcf_per_branch(tree = tree, vcf = vcf, ref_genome = ref_genome)
branch_grl = convert_vcf_to_granges(branch_vcf_list = branch_vcf, ref_genome = ref_genome)
branch_mm = mut_matrix(branch_grl, ref_genome)

# Preprocess bulk samples

vcf_bulk_snvs_LN <- vcf_bulk_LN[isSNV(vcf_bulk_LN)] # only SNVs
autosomes <- as.character(1:22) 
vcf_bulk_autosomal_LN <- vcf_bulk_snvs_LN[seqnames(rowRanges(vcf_bulk_snvs_LN)) %in% autosomes] # only autosomes 
vaf_matrix_LN <- geno(vcf_bulk_autosomal_LN)$VAF
vaf_matrix_LN <- vaf_matrix_LN[, 2, drop = FALSE]  # keep as a matrix

# Function to extract variant key

variant_key <- function(vcf) {
  rr <- rowRanges(vcf)
  paste0(
    seqnames(rr), ":", start(rr), "_",
    as.character(mcols(rr)$REF), ">",
    as.character(unlist(mcols(rr)$ALT))
  )
}

# Get variant keys and VAFs for the bulk

bulk_keys <- variant_key(vcf_bulk_autosomal_LN)
bulk_vaf_matrix <- geno(vcf_bulk_autosomal_LN)$VAF
bulk_vafs <- as.numeric(bulk_vaf_matrix[, 2]) 

# Create a named vector of VAFs for lookup

names(bulk_vafs) <- bulk_keys

# Prepare result list

shared_stats <- lapply(names(branch_vcf), function(branch_name) {
  branch <- branch_vcf[[branch_name]]
  branch_keys <- variant_key(branch)
  
  # Get intersection
  shared <- intersect(branch_keys, names(bulk_vafs))
  
  # Get VAFs of shared variants from bulk
  shared_vafs <- bulk_vafs[shared]
  
  # Return as a list or data frame row
  list(
    branch = branch_name,
    n_shared = length(shared),
    shared_vafs = shared_vafs
  )
})

# Convert to data frame with summary info

shared_summary_df <- do.call(rbind, lapply(shared_stats, function(x) {
  data.frame(branch = x$branch, n_shared = x$n_shared)
}))

# Optionally view full VAF distributions per branch

names(shared_stats) <- sapply(shared_stats, function(x) x$branch)

# Example: plot VAF distributions per group

group1 <- c("h")
group2 <- c("p", "q", "a")
group3 <- c("G", "F", "E", "C", "B", "A", "v", "u", "x", "t", "M", "L", "O", "Q", "S", "V", "U", "Z", "Y")
group4 <- c("o", "m")
group5 <- c("j", "i", "l", "n")
group6 <- c("g", "b", "c")

# Helper to combine VAFs from a list of branch names
get_vafs_for_group <- function(branch_names, stats_list) {
  unlist(lapply(branch_names, function(branch) {
    if (!is.null(stats_list[[branch]])) {
      stats_list[[branch]]$shared_vafs
    } else {
      numeric(0)  # In case branch is missing
    }
  }))
}

vafs_group1 <- get_vafs_for_group(group1, shared_stats)
vafs_group2 <- get_vafs_for_group(group2, shared_stats)
vafs_group3 <- get_vafs_for_group(group3, shared_stats)
vafs_group4 <- get_vafs_for_group(group4, shared_stats)
vafs_group5 <- get_vafs_for_group(group5, shared_stats)
vafs_group6 <- get_vafs_for_group(group6, shared_stats)

# Build one long data frame for ggplot
plot_df <- rbind(
  data.frame(VAF = vafs_group1, Group = "BL-Trunk"),
  data.frame(VAF = vafs_group2, Group = "BL-Intermediate-LN"),
  data.frame(VAF = vafs_group3, Group = "BL-Private-LN"),
  #data.frame(VAF = vafs_group4, Group = "BL-Intermediate-BM"), # no mutations
  data.frame(VAF = vafs_group5, Group = "BL-Private-BM"),
  data.frame(VAF = vafs_group6, Group = "WT")
)

tree_mut_counts <- data.frame(
  branch = tree@data$branch_id,
  total_muts = tree@data$branch_length  
)

get_shared_vs_total_stats <- function(branches, stats_list, tree_df) {
  shared <- sum(sapply(branches, function(b) length(stats_list[[b]]$shared_vafs)))
  total <- sum(tree_df$total_muts[match(branches, tree_df$branch)], na.rm = TRUE)
  
  data.frame(
    BranchGroup = paste(branches, collapse = "-"),  # renamed to avoid conflict
    Shared = shared,
    Total = total,
    Label = paste0(shared, "/", total, " = ", round(100 * shared / total, 1), "%")
  )
}

summary_stats <- rbind(
  cbind(get_shared_vs_total_stats(group1, shared_stats, tree_mut_counts), Group = "BL-Trunk"),
  cbind(get_shared_vs_total_stats(group2, shared_stats, tree_mut_counts), Group = "BL-Intermediate-LN"),
  cbind(get_shared_vs_total_stats(group3, shared_stats, tree_mut_counts), Group = "BL-Private-LN"),
  cbind(get_shared_vs_total_stats(group4, shared_stats, tree_mut_counts), Group = "BL-Intermediate-BM"),
  cbind(get_shared_vs_total_stats(group5, shared_stats, tree_mut_counts), Group = "BL-Private-BM"),
  cbind(get_shared_vs_total_stats(group6, shared_stats, tree_mut_counts), Group = "WT")
)

plot_df$Group <- factor(plot_df$Group, levels = c("BL-Trunk", "BL-Intermediate-LN", "BL-Private-LN", "BL-Intermediate-BM","BL-Private-BM", "WT"))
summary_stats$Group <- factor(summary_stats$Group, levels = c("BL-Trunk", "BL-Intermediate-LN", "BL-Private-LN", "BL-Intermediate-BM","BL-Private-BM", "WT"))

group_cols <- c(
  "BL-Trunk" = "#1D3557",
  "BL-Intermediate-LN" = "#0A9086",
  "BL-Private-LN" = "#A62639",
  "BL-Intermediate-BM" = "#0A9086",
  "BL-Private-BM" = "#A62639",
  "WT" = "#D2BD96"
)


p4a <- ggplot(plot_df, aes(x = VAF)) +
  geom_histogram(aes(fill = Group), bins = 10, color = "white") +
  facet_wrap(~Group, scales = "free_y") +
  scale_fill_manual(values = group_cols) +
  coord_cartesian(xlim = c(0, 1)) +
  geom_text(
    data = summary_stats,
    aes(x = 0.8, y = Inf, label = Label),
    inherit.aes = FALSE,
    vjust = 2, size = 2.5, fontface = "bold"
  ) +
  theme_classic() +
  guides(fill = "none") +  # legend not needed since facets are labeled
  labs(title = "Bulk VAF Distributions of Shared Mutations",
       x = "Bulk VAF", y = "Count")

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_bulk_VAF_distribution_histogram_LN.pdf",
  plot = p4a,
  width = 8,
  height = 6
)

p5a <- ggplot(plot_df, aes(y = VAF)) +
  geom_violin(aes(x = 1, fill = Group), color = "white") +
  facet_wrap(~ Group, scales = "free_y") +
  scale_fill_manual(values = group_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_text(
    data = summary_stats,
    aes(x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,label = Label),
    inherit.aes = FALSE,
    vjust = 2, size = 3, fontface = "bold"
  ) +
  scale_x_continuous(breaks = NULL) +
  theme_classic() +
  guides(fill = "none") +
  labs(title = "Bulk VAF Distributions of Shared Mutations",
       x = NULL, y = "Bulk VAF")

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_bulk_VAF_distribution_violin_LN.pdf",
  plot = p5a,
  width = 8,
  height = 6
)

# Look within intermediate branches to see if the clade size matches median VAF 

intermediate_list <- shared_stats[intersect(names(shared_stats), group2)]

# Stack into a long data.frame, skipping branches with no shared_vafs
dfs <- lapply(names(intermediate_list), function(b) {
  vals <- intermediate_list[[b]]$shared_vafs
  vals <- vals[!is.na(vals)]              # drop NAs
  if (length(vals) == 0) return(NULL)     # skip empties
  data.frame(
    Branch = rep(b, length(vals)),
    VAF    = vals,
    stringsAsFactors = FALSE
  )
})

intermediate_df <- do.call(rbind, dfs)

# In case all were empty:
if (is.null(intermediate_df) || nrow(intermediate_df) == 0) {
  stop("No shared VAFs found for the Intermediate branches in shared_stats.")
}

p6a <- ggplot(intermediate_df %>% group_by(Branch) %>% mutate(nuniq = n_distinct(VAF)) %>% ungroup(), aes(x = 1, y = VAF)) +
  geom_violin(data = ~ subset(.x, nuniq >= 2), fill = "#0A9086", color = "white", trim = FALSE) +
  geom_point(data = ~ subset(.x, nuniq < 2), position = position_jitter(width = 0.04, height = 0), alpha = 0.6, size = 1) +
  facet_wrap(~ Branch, scales = "free_y") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = NULL) +
  theme_classic() +
  labs(title = "Intermediate branches: VAF distributions of shared mutations", x = NULL, y = "Bulk VAF")

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_bulk_VAF_distribution_intermediate_branches_violin_LN.pdf",
  plot = p6a,
  width = 8,
  height = 6
)

# now repeat for BM bulk

tumor_sample <- "PMABM000FZO"
bulk_df_BM <- data.frame()

if (tumor_sample %in% samples(header(vcf_bulk_BM_snv_autosomal))) {
  ad <- sapply(geno(vcf_bulk_BM_snv_autosomal)$AD[, tumor_sample], "[[", 2)  # Alt depth
  dp <- geno(vcf_bulk_BM_snv_autosomal)$DP[, tumor_sample]  # Total depth
  vaf <- ad / dp
  
  bulk_df_BM <- data.frame(variant = names(vaf), vaf = as.numeric(vaf), bulk = tumor_sample)
}
bulk_df_BM_filtered <- bulk_df_BM %>%
  filter(vaf >= 0.15)

bulk_name <- unique(bulk_df_BM_filtered$bulk)[1]

vaf_BM <- bulk_df_BM_filtered |>
  dplyr::select(variant, vaf) |>
  tibble::column_to_rownames("variant") |>
  setNames(bulk_name) |>
  as.matrix()

rownames(vaf_BM) <- sub("^chr", "", rownames(vaf_BM))
vaf_BM <- setNames(as.numeric(vaf_BM[, 1]), rownames(vaf_BM))
names(vaf_BM) <- gsub("/", ">", names(vaf_BM), fixed = TRUE)

# Prepare result list

shared_stats <- lapply(names(branch_vcf), function(branch_name) {
  branch <- branch_vcf[[branch_name]]
  branch_keys <- variant_key(branch)
  
  # Get intersection
  shared <- intersect(branch_keys, names(vaf_BM))
  
  # Get VAFs of shared variants from bulk
  shared_vafs <- vaf_BM[shared]
  
  # Return as a list or data frame row
  list(
    branch = branch_name,
    n_shared = length(shared),
    shared_vafs = shared_vafs
  )
})

# Convert to data frame with summary info

shared_summary_df <- do.call(rbind, lapply(shared_stats, function(x) {
  data.frame(branch = x$branch, n_shared = x$n_shared)
}))

# Optionally view full VAF distributions per branch

names(shared_stats) <- sapply(shared_stats, function(x) x$branch)

# Helper to combine VAFs from a list of branch names
get_vafs_for_group <- function(branch_names, stats_list) {
  unlist(lapply(branch_names, function(branch) {
    if (!is.null(stats_list[[branch]])) {
      stats_list[[branch]]$shared_vafs
    } else {
      numeric(0)  # In case branch is missing
    }
  }))
}

vafs_group1 <- get_vafs_for_group(group1, shared_stats)
vafs_group2 <- get_vafs_for_group(group2, shared_stats)
vafs_group3 <- get_vafs_for_group(group3, shared_stats)
vafs_group4 <- get_vafs_for_group(group4, shared_stats)
vafs_group5 <- get_vafs_for_group(group5, shared_stats)
vafs_group6 <- get_vafs_for_group(group6, shared_stats)

# Build one long data frame for ggplot
plot_df <- rbind(
  data.frame(VAF = vafs_group1, Group = "BL-Trunk"),
  data.frame(VAF = vafs_group2, Group = "BL-Intermediate-LN"),
  #data.frame(VAF = vafs_group3, Group = "BL-Private-LN"), # no mutations
  data.frame(VAF = vafs_group4, Group = "BL-Intermediate-BM"), 
  data.frame(VAF = vafs_group5, Group = "BL-Private-BM")
  #data.frame(VAF = vafs_group6, Group = "WT") # no mutations
)

tree_mut_counts <- data.frame(
  branch = tree@data$branch_id,
  total_muts = tree@data$branch_length  
)

get_shared_vs_total_stats <- function(branches, stats_list, tree_df) {
  shared <- sum(sapply(branches, function(b) length(stats_list[[b]]$shared_vafs)))
  total <- sum(tree_df$total_muts[match(branches, tree_df$branch)], na.rm = TRUE)
  
  data.frame(
    BranchGroup = paste(branches, collapse = "-"),  # renamed to avoid conflict
    Shared = shared,
    Total = total,
    Label = paste0(shared, "/", total, " = ", round(100 * shared / total, 1), "%")
  )
}

summary_stats <- rbind(
  cbind(get_shared_vs_total_stats(group1, shared_stats, tree_mut_counts), Group = "BL-Trunk"),
  cbind(get_shared_vs_total_stats(group2, shared_stats, tree_mut_counts), Group = "BL-Intermediate-LN"),
  cbind(get_shared_vs_total_stats(group3, shared_stats, tree_mut_counts), Group = "BL-Private-LN"),
  cbind(get_shared_vs_total_stats(group4, shared_stats, tree_mut_counts), Group = "BL-Intermediate-BM"),
  cbind(get_shared_vs_total_stats(group5, shared_stats, tree_mut_counts), Group = "BL-Private-BM"),
  cbind(get_shared_vs_total_stats(group6, shared_stats, tree_mut_counts), Group = "WT")
)

plot_df$Group <- factor(plot_df$Group, levels = c("BL-Trunk", "BL-Intermediate-LN", "BL-Private-LN", "BL-Intermediate-BM","BL-Private-BM", "WT"))
summary_stats$Group <- factor(summary_stats$Group, levels = c("BL-Trunk", "BL-Intermediate-LN", "BL-Private-LN", "BL-Intermediate-BM","BL-Private-BM", "WT"))

group_cols <- c(
  "BL-Trunk" = "#1D3557",
  "BL-Intermediate-LN" = "#0A9086",
  "BL-Private-LN" = "#A62639",
  "BL-Intermediate-BM" = "#0A9086",
  "BL-Private-BM" = "#A62639",
  "WT" = "#D2BD96"
)


p4b <- ggplot(plot_df, aes(x = VAF)) +
  geom_histogram(aes(fill = Group), bins = 10, color = "white") +
  facet_wrap(~Group, scales = "free_y") +
  scale_fill_manual(values = group_cols) +
  coord_cartesian(xlim = c(0, 1)) +
  geom_text(
    data = summary_stats,
    aes(x = 0.8, y = Inf, label = Label),
    inherit.aes = FALSE,
    vjust = 2, size = 2.5, fontface = "bold"
  ) +
  theme_classic() +
  guides(fill = "none") +  # legend not needed since facets are labeled
  labs(title = "Bulk VAF Distributions of Shared Mutations",
       x = "Bulk VAF", y = "Count")

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_bulk_VAF_distribution_histogram_BM.pdf",
  plot = p4b,
  width = 8,
  height = 6
)

p5b <- ggplot(plot_df, aes(y = VAF)) +
  geom_violin(aes(x = 1, fill = Group), color = "white") +
  facet_wrap(~ Group, scales = "free_y") +
  scale_fill_manual(values = group_cols) +
  coord_cartesian(ylim = c(0, 1)) +
  geom_text(
    data = summary_stats,
    aes(x = Inf, y = Inf, hjust = 1.05, vjust = 1.4,label = Label),
    inherit.aes = FALSE,
    vjust = 2, size = 3, fontface = "bold"
  ) +
  scale_x_continuous(breaks = NULL) +
  theme_classic() +
  guides(fill = "none") +
  labs(title = "Bulk VAF Distributions of Shared Mutations",
       x = NULL, y = "Bulk VAF")

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_bulk_VAF_distribution_violin_BM.pdf",
  plot = p5b,
  width = 8,
  height = 6
)

# Look within intermediate branches to see if the clade size matches median VAF 

intermediate_list <- shared_stats[intersect(names(shared_stats), group4)]

# Stack into a long data.frame, skipping branches with no shared_vafs
dfs <- lapply(names(intermediate_list), function(b) {
  vals <- intermediate_list[[b]]$shared_vafs
  vals <- vals[!is.na(vals)]              # drop NAs
  if (length(vals) == 0) return(NULL)     # skip empties
  data.frame(
    Branch = rep(b, length(vals)),
    VAF    = vals,
    stringsAsFactors = FALSE
  )
})

intermediate_df <- do.call(rbind, dfs)

# In case all were empty:
if (is.null(intermediate_df) || nrow(intermediate_df) == 0) {
  stop("No shared VAFs found for the Intermediate branches in shared_stats.")
}

p6b <- ggplot(intermediate_df %>% group_by(Branch) %>% mutate(nuniq = n_distinct(VAF)) %>% ungroup(), aes(x = 1, y = VAF)) +
  geom_violin(data = ~ subset(.x, nuniq >= 2), fill = "#0A9086", color = "white", trim = FALSE) +
  geom_point(data = ~ subset(.x, nuniq < 2), position = position_jitter(width = 0.04, height = 0), alpha = 0.6, size = 1) +
  facet_wrap(~ Branch, scales = "free_y") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = NULL) +
  theme_classic() +
  labs(title = "Intermediate branches: VAF distributions of shared mutations", x = NULL, y = "Bulk VAF")

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_bulk_VAF_distribution_intermediate_branches_violin_BM.pdf",
  plot = p6b,
  width = 8,
  height = 6
)


# Plot signatures on tree

# Get signatures 
all_signatures = get_known_signatures()
pta_v1_sig = read.table("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_lymphoma/1_Input/WGS/PTA_Artefact_Signature.txt", sep = "\t", header = T)
pta_v1_sig = as.matrix(pta_v1_sig)
PTA_v1 <- as.numeric(pta_v1_sig[,"PTA"])
PTA_v1 <- PTA_v1[!is.na(PTA_v1)]
pta_v2_sig = read.table("~/hpc/pmc_vanboxtel/resources/signatures/PTAv2_Artefact_Signature.txt", sep = "\t", header = T) 
pta_v2_sig = as.matrix(pta_v2_sig) 
PTA_v2 <- as.numeric(pta_v2_sig[,"PTAv2"])
sbsblood <- read.table("~/Downloads/sigfit_cosmic3_bloodsig_Aug2020.txt", sep = "\t", header = T)
sbsblood = as.matrix(sbsblood)
SBSblood <- as.numeric(sbsblood[,"Signature.Blood"])
signatures <- cbind(SBSblood, all_signatures, PTA_v1, PTA_v2)

# Refitting part without PTA artefact signature 

sub_sig <- signatures[, c("SBS1", "SBS9", "SBS17b", "SBS18", "SBSblood")] # from de novo extraction

contribution <- fit_to_signatures_strict_tree(mut_matrix = branch_mm, signatures = sub_sig, max_delta = 0.002, remove_min = 0)

pal <- c("#D2BD96", "#0A9086", "#B3B3B3", "#A62639", "#1D3557", "#6C5B7B")
names(pal) <- c("SBS1","SBS9","SBS17b","SBS18","SBSblood","PTA_v1")
sig_cols <- pal[rownames(contribution)]

plot_contribution(contribution, palette = pal)

# check if all branches are explained well with these signatures: cosine < 0.85 with > 200 mutations would suggest you miss a mutation
p7 <- check_reconstructed_cosine(contribution, branch_mm, sub_sig, tree) +
  ylim(0,1)

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/PRN4_reconstructed_cosine_without_PTA.pdf",
  plot = p7,
  width = 8,
  height = 6
)

# missing some signture(s) so will include PTA artifact signature

sub_sig <- signatures[, c("SBS1", "SBS9", "SBS17b", "SBS18", "SBSblood", "PTA_v1")] 

contribution <- fit_to_signatures_strict_tree(mut_matrix = branch_mm, signatures = sub_sig, max_delta = 0.002, remove_min = 0)

pal <- c("#D2BD96", "#0A9086", "#B3B3B3", "#A62639", "#1D3557", "#6C5B7B")
names(pal) <- c("SBS1","SBS9","SBS17b","SBS18","SBSblood","PTA_v1")
sig_cols <- pal[rownames(contribution)]

plot_contribution(contribution, palette = pal)

# check if all branches are explained well with these signatures: cosine < 0.85 with > 200 mutations would suggest you miss a mutation
p8 <- check_reconstructed_cosine(contribution, branch_mm, sub_sig, tree) +
  ylim(0,1)

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/PRN4_reconstructed_cosine_with_PTA.pdf",
  plot = p8,
  width = 8,
  height = 6
)

# most branches >200 mutations are above the 0.85 threshold so keep PTA signature

# add signature contributions to your tree
tree = add_contribution(tree, contribution = contribution) # if you already did signature_fitting

p9 <- plot_tree_contribution_bars_new_alex(
  tree = tree,
  signatures = sub_sig,
  mut_matrix = branch_mm,
  signature_colors = sig_cols,
  title = "PRN4",
  bar_height = 0.02,
  scaling = 0.8,
  branch_color = "#000000",
  x_limits = c(0, 2600),
  x_breaks = seq(0, 2600, by = 200)
)

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/PRN4_signatures_in_tree_with_PTA.pdf",
  plot = p9,
  width = 8,
  height = 6
)

# Convert the matrix to a data frame first
contribution_df <- as.data.frame(contribution)

# Get the column names of the data frame
column_names <- colnames(contribution_df)

# Assign groups based on column names
group_assignment <- case_when(
  column_names %in% group1 ~ "BL-Trunk",
  column_names %in% group2 ~ "BL-Intermediate-LN",
  column_names %in% group3 ~ "BL-Private-LN",
  column_names %in% group4 ~ "BL-Intermediate-BM",
  column_names %in% group5 ~ "BL-Private-BM",
  column_names %in% group6 ~ "WT"
)

# Add the group information to the data frame
contribution_with_group <- data.frame(t(contribution_df), group = group_assignment)
contribution_with_group <- contribution_with_group %>%
  pivot_longer(cols = -group, names_to = "sample", values_to = "value")

# Now, group by 'group' and summarise the contributions
contribution_grouped <- contribution_with_group %>%
  group_by(group, sample) %>%
  summarise(total_contribution = sum(value, na.rm = TRUE), .groups = 'drop')

# Pivot wider to have groups as columns and signatures (samples) as rows
contribution_wide <- contribution_grouped %>%
  pivot_wider(names_from = group, values_from = total_contribution, values_fill = list(total_contribution = 0))

# Add BL
contribution_wide <- contribution_wide |>
  mutate(
    BL = coalesce(`BL-Intermediate-LN`, 0) + coalesce(`BL-Intermediate-BM`, 0) + coalesce(`BL-Private-LN`, 0) + coalesce(`BL-Private-BM`, 0)
  )

# Re-order

contribution_wide <- contribution_wide |>
  dplyr::select(sample, WT, BL, "BL-Trunk", "BL-Intermediate-LN", "BL-Private-LN", "BL-Intermediate-BM", "BL-Private-BM" )

# Convert the result into a matrix, using sample names as row names
contribution_matrix <- as.matrix(contribution_wide %>% column_to_rownames("sample"))

# View the resulting matrix
contribution_matrix

p10 <- plot_contribution(contribution_matrix[,], palette = sig_cols) # All signatures including PTA

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/PRN4_signatures_per_group_with_PTA.pdf",
  plot = p10,
  width = 8,
  height = 6
)

p11 <- plot_contribution(contribution_matrix[2:6,], palette = sig_cols) # without PTA

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/PRN4_signatures_per_group_without_PTA.pdf",
  plot = p11,
  width = 8,
  height = 6
)

p11_sub <- plot_contribution(contribution_matrix[2:6,3:7], palette = sig_cols) # without PTA

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/PRN4_signatures_per_group_without_PTA_sub.pdf",
  plot = p11_sub,
  width = 8,
  height = 6
)

# Convert to long format

df_long <- melt(contribution_matrix[2:6,], varnames = c("Signature", "Sample"), value.name = "Contribution") %>%
  group_by(Sample) %>%
  mutate(Prop = Contribution / sum(Contribution)) %>%
  ungroup()

p11b <- ggplot(df_long, aes(x = Sample, y = Prop, fill = Signature)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = sig_cols) +
  scale_y_continuous(limits = c(0, 1)) +     # y-axis from 0 to 1
  labs(title = "Signature contribution per sample (no PTA)",
       x = "",
       y = "")

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/PRN4_signatures_per_group_without_PTA_separate_bars.pdf",
  plot = p11b,
  width = 8,
  height = 6
)

p11b_sub <- ggplot(transform(subset(df_long, Sample %in% c("BL-Trunk","BL-Intermediate-LN","BL-Private-LN", "BL-Intermediate-BM","BL-Private-BM")), Sample = factor(Sample, levels = c("BL-Trunk","BL-Intermediate-LN","BL-Private-LN", "BL-Intermediate-BM","BL-Private-BM"))), aes(x = Sample, y = Prop, fill = Signature)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = sig_cols) +
  scale_y_continuous(limits = c(0, 1)) +     # y-axis from 0 to 1
  labs(title = "Signature contribution per sample (no PTA)",
       x = "",
       y = "")

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/PRN4_signatures_per_group_without_PTA_separate_bars_sub.pdf",
  plot = p11b_sub,
  width = 8,
  height = 6
)

# Get relative proportion of C>TpG sites

# To do that first need to correct branch_mm for PTAsigContribution

pta_contri_df_wide <- as.data.frame(t(contribution_df["PTA_v1", , drop = FALSE]))
colnames(pta_contri_df_wide) <- "PTA_v1"

subtract_ptasig <- function(input_mutmat) {
  
  # get the absolute PTA contributions
  pta_contris <- pta_contri_df_wide[colnames(input_mutmat),]
  # reconstruct mutmat of PTA contribution
  reconstructed_PTAmutmats <- as.matrix(PTA_v1) %*% as.matrix(t(pta_contris))
  # subtract
  output_mutmat <- input_mutmat - reconstructed_PTAmutmats
  output_mutmat[output_mutmat < 0] <- 0 
  
  return(output_mutmat)
}

branch_mm_corrected <- subtract_ptasig(branch_mm)

sample_groups <- c(
  setNames(rep("BL-Trunk", length(group1)), group1),
  setNames(rep("BL-Intermediate-LN", length(group2)), group2),
  setNames(rep("BL-Private-LN", length(group3)), group3),
  setNames(rep("BL-Intermediate-BM", length(group4)), group4),
  setNames(rep("BL-Intermediate-BM", length(group5)), group5),
  setNames(rep("WT", length(group6)), group6)
)

# Add merged BL group
BL_samples <- c(group2, group3, group4, group5)
sample_groups <- c(sample_groups, setNames(rep("BL", length(BL_samples)), BL_samples))

# Subset matrix to only relevant samples

branch_sel <- branch_mm_corrected[, names(sample_groups), drop = FALSE]

# Sum counts per group

branch_by_group <- rowsum(t(branch_sel), group = sample_groups)  # groups in rows
branch_by_group <- t(branch_by_group)                            # groups as columns

# Find CpG C>T channels

ix_cpg_ct <- grepl("\\[C>T\\]G$", rownames(branch_by_group))

# Compute proportions

cpg_ct <- colSums(branch_by_group[ix_cpg_ct, , drop = FALSE])
totals <- colSums(branch_by_group)
prop <- cpg_ct / totals

# Prepare data frame for plotting

df <- data.frame(Group = names(prop), Proportion = as.numeric(prop))

# Plot

p11c <- ggplot(df, aes(x = Group, y = Proportion)) +
  geom_col(fill = "#D2BD96") +
  labs(
    title = "Fraction of C>T at CpG sites per group",
    x = "",
    y = ""
  ) +
  ylim(0, 1)

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/PRN4_CTpG_fraction_custom_groups.pdf",
  plot = p11c,
  width = 6,
  height = 4.5
)

# remove PTA signature contribution and re-plot absolute tree 

ggtree(tree@phylo)
metadata <- tree@data
print(metadata)
x <- as_tibble(tree)
phylo_tree <-as.phylo(x)
x[28,] #this row contains NAs (node 15-to-15 connection)
x<- x[-28,] #deleting NA-filled row
phylo_tree$edge.length<-unname(x$branch_length)
phylo_tree$edge
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label,  # Must match tree$tip.label
  age = rep(4.1,52)  # Replace these with your actual ages
)
phylo_tree$agedf<-agedf
ggtree(phylo_tree)

# read csv with signature per branch info
sig_contri <- contribution
sbs_minus_pta <- colSums(sig_contri[c("SBSblood", "SBS1", "SBS9", "SBS18", "SBS17b"), ])
names(sbs_minus_pta) <- colnames(sig_contri)
missing_branches <- setdiff(x$branch_id, names(sbs_minus_pta))
sbs_minus_pta[missing_branches] <- 0
x$branch_length_minus_PTA <- sbs_minus_pta[ x$branch_id ]
x$branch_length_minus_PTA[is.na(x$branch_length_minus_PTA)] <- 0
x$branch_length_minus_PTA <- as.integer(x$branch_length_minus_PTA)

# Define the product sensitivity calculator
calculate_product_sensitivity <- function(samples_str, sensitivity_df) {
  sample_list <- strsplit(samples_str, "\\|")[[1]]
  sensitivities <- sensitivity_df$Callable_fraction[sensitivity_df$Sample_name %in% sample_list]
  if (length(sensitivities) == 0) return(1)  # fallback to avoid division by zero
  return(1 - prod(1 - sensitivities, na.rm = TRUE))
}

x <- x %>%
  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_PRN4_filtered),
         corr_branch_lengths = x$branch_length_minus_PTA / product_sensitivity)

x <- x %>%
  mutate(corr_branch_lengths = round(corr_branch_lengths))

phylo_tree$edge.length<-unname(x$branch_length_minus_PTA)
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label,  # Must match tree$tip.label
  age = rep(4.1,52)  # Replace these with your actual ages
)
phylo_tree$agedf<-agedf
ggtree(phylo_tree, )

tree@data$branch_length <- x$corr_branch_lengths[ match(tree@data$branch_id, x$branch_id) ]


# add signature contributions to your tree
phylo_tree = add_contribution(tree, contribution = contribution[1:5,]) # if you already did signature_fitting

p12 <- plot_tree_contribution_bars_new_alex(
  tree = tree, contribution = contribution[1:5,], signature_colors = sig_cols,
  title = "PRN4", bar_height = 0.02, scaling = 0.8, add_branch_length = F,
  branch_color = "#000000",
  x_limits = c(0, 2200),
  x_breaks = seq(0, 2200, by = 200)
)

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/PRN4_signatures_in_tree_without_PTA_callable_loci_corrected.pdf",
  plot = p12,
  width = 8,
  height = 6
)

# Plot final corrected tree
# Merge the tree tip data with your annotation data
tip_data <- input_df_PRN4 %>%
  mutate(
    color_group = ifelse(Myc_translocation_IGV == "Yes", "#4378bd", "#e7872b"),
    shape_group = ifelse(Biopsy_type == "LN", 16, 17)  # 16 = filled circle, 17 = filled triangle
  )

# Plot tree with shapes and colors
p13 <- plot_gg_tree(
  tree,
  add_branch_length = FALSE,
  add_bootstrap = FALSE,
  add_tip_label = FALSE,
  add_title = ""
) %<+% tip_data +   # attach annotation data to tree
  
  geom_tree(color = "grey") +
  
  # Tip points colored by Myc_translocation_IGV and shaped by Biopsy_type
  geom_tippoint(aes(color = Myc_translocation_IGV, shape = Biopsy_type), size = 3) +
  
  # Colors for Myc_translocation_IGV
  scale_color_manual(values = c("Yes" = "#4378bd", "No" = "#e7872b")) +
  
  # Shapes for Biopsy_type
  scale_shape_manual(values = c("LN" = 16, "BM" = 17)) +
  
  scale_x_continuous(
    name = "Mutation burden (SNVs)",
    limits = c(0, 2200),
    breaks = seq(0, 2200, by = 200)
  ) +
  
  theme(
    axis.text.x = element_text(size = 10, angle = 90,
                               vjust = 0.5, hjust = 1, margin = margin(t = 6)),
    axis.title.x = element_text(size = 12, angle = 180,
                                vjust = 0.5, margin = margin(t = 6))
  )

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/PRN4_absolute_tree_without_PTA_corrected_for_callable_loci.pdf",
  plot = p13,
  width = 8,
  height = 6
)

