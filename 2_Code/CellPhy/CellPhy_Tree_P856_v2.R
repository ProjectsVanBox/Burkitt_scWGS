################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at CellPhy Tree, branch signatures and Ultrametric tree for Donor P856
# Author: Alexander Steemers
# Date: August 2025
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
source("~/surfdrive/Shared/pmc_vanboxtel/personal/asteemers/R_functions/check_reconstructed_cosine_Alex.R")
source("~/surfdrive/Shared/pmc_vanboxtel/personal/asteemers/R_functions/plot_gg_tree_alex.R")
source("~/surfdrive/Shared/pmc_vanboxtel/personal/asteemers/R_functions/correct_branches.R")
source("~/surfdrive/Shared/pmc_vanboxtel/personal/asteemers/R_functions/calculate_product_sensitivity.R")
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')
theme_set(theme_classic())
source("~/surfdrive/Shared/pmc_vanboxtel/general/2_Bioinformatics/colors/Jurrians_colors.R")

# Load tree object and vcf files

tree = readRDS("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/P856/Filtered_samples_v5/CPW_01/TreeObject0.1.RDS")
vcf = VariantAnnotation::readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/P856/Filtered_samples_v5/P856.vep.sub.SNV.autosomal.noBulksNoBlacklist.vcf")
vcf_bulk = readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P856/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P856_bulk.vep.SMuRF.filtered.sorted.vcf.gz")

# Prepare tree

tree = prepare_tree(tree)

# Plot bare tree

p1 <- plot_gg_tree_base(tree, add_tip_label = T, add_title = "")

# Save as PDF

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_absolute_tree_with_branch_info.pdf",
  plot = p1,
  width = 8,
  height = 6
)

p2 <- plot_gg_tree(tree, add_branch_length = F, add_bootstrap = F, add_tip_label = F,add_title = "") +
  geom_tree(color = "grey") +  # Set branch color to grey
  scale_x_continuous(name = "Mutation burden (SNVs)",
                     limits = c(0, 3400),
                     breaks = seq(0, 3400, by = 200)) +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 6)),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5, margin = margin(t = 6) )
  )

# Save as PDF

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_absolute_tree.pdf",
  plot = p2,
  width = 8,
  height = 6
)

p2b <- plot_gg_tree_alex(tree ,add_branch_length = F, add_bootstrap = T, add_tip_label = F,add_title = "", add_x = FALSE) +
  geom_tree(color = "grey") +  # Set branch color to grey
  scale_x_continuous(name = "Mutation burden (SNVs)",
                     limits = c(0, 3400),
                     breaks = seq(0, 3400, by = 200)) +
  theme_tree2() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 6)),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5, margin = margin(t = 6) )
  )

# Save as PDF

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_absolute_tree_with_bootstraps.pdf",
  plot = p2b,
  width = 8,
  height = 6
)

# Import callable loci for patient

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') #dataframe
input_df_P856 <- input_df[input_df$Novogene_ID == "P856" & input_df$ResolveDNA_version %in% c("v1", "v2.0", "v2"),]
below_curve_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv")
low_call_frac_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv")
fail_vaf_df <-  read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt")
filtered_samples <- unique(c(low_call_frac_df$Sample_name, below_curve_df$Sample_name, fail_vaf_df$samplename))  # samples that didn't pass QC
input_df_P856_filtered <- input_df_P856 %>% filter(!Sample_name %in% filtered_samples)
input_df_P856_filtered$Callable_Loci <- as.numeric(input_df_P856_filtered$Callable_Loci)

# Correct the branches

corrected_tree <- correct_branches(tree, input_df_P856_filtered)
corrected_tree@data$branch_length <- round(corrected_tree@data$branch_length)

# Plot corrected tree

tip_data <- input_df_P856 %>%
  mutate(color_group = ifelse(Myc_translocation_IGV == "Yes", "#4378bd", "#e7872b"))

p3 <- plot_gg_tree(corrected_tree, add_branch_length = F, add_bootstrap = F, add_tip_label = F,add_title = "") %<+% tip_data + 
  geom_tree(color = "grey") +  # Set branch color to grey
  geom_tippoint(aes(color = Myc_translocation_IGV), size = 1) +  # map color to Yes/No
  scale_color_manual(values = c("Yes" = "#4378bd", "No" = "#e7872b")) +
  scale_x_continuous(name = "Mutation burden (SNVs)",
                     limits = c(0, 3600),
                     breaks = seq(0, 3600, by = 200)) +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 6)),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5, margin = margin(t = 6) )
  )

# Save as PDF

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_absolute_tree_corrected_for_callable_loci.pdf",
  plot = p3,
  width = 8,
  height = 6
)

p3b <- plot_gg_tree_alex(corrected_tree ,add_branch_length = F, add_bootstrap = T, add_tip_label = F,add_title = "", add_x = FALSE) +
  geom_tree(color = "grey") +  # Set branch color to grey
  scale_x_continuous(name = "Mutation burden (SNVs)",
                     limits = c(0, 3600),
                     breaks = seq(0, 3600, by = 200)) +
  theme_tree2() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1, margin = margin(t = 6)),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5, margin = margin(t = 6) )
  )

# Save as PDF

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_absolute_tree_with_bootstraps_corrected_for_callable_loci.pdf",
  plot = p3b,
  width = 8,
  height = 6
)

# VCF per branch

branch_vcf = extract_vcf_per_branch(tree = tree, vcf = vcf, ref_genome = ref_genome)
saveRDS(branch_vcf, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_vcf_P856.rds")
branch_grl = convert_vcf_to_granges(branch_vcf_list = branch_vcf, ref_genome = ref_genome)
saveRDS(branch_grl, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_grl_P856.rds")
branch_mm = mut_matrix(branch_grl, ref_genome)

# Preprocess bulk samples

vcf_bulk_snvs <- vcf_bulk[isSNV(vcf_bulk)] # only SNVs
autosomes <- as.character(1:22) 
vcf_bulk_autosomal <- vcf_bulk_snvs[seqnames(rowRanges(vcf_bulk_snvs)) %in% autosomes] # only autosomes 
vaf_matrix <- geno(vcf_bulk_autosomal)$VAF
vaf_matrix_BM <- vaf_matrix[, 1, drop = FALSE]  # keep as a matrix
vaf_matrix_PL <- vaf_matrix[, 3, drop = FALSE]  # keep as a matrix

# Function to extract variant key

variant_key <- function(vcf) {
  rr <- rowRanges(vcf)
  paste0(
    seqnames(rr), ":", start(rr), "_",
    as.character(mcols(rr)$REF), ">",
    as.character(unlist(mcols(rr)$ALT))
  )
}

# Get variant keys and VAFs for the Pleura bulk 

bulk_keys <- variant_key(vcf_bulk_autosomal)
vaf_mat   <- geno(vcf_bulk_autosomal)$VAF
vaf_PL      <- as.numeric(vaf_mat[, 3])
mask      <- !is.na(vaf_PL) & vaf_PL >= 0.15

bulk_vafs_PL <- vaf_PL[mask]
names(bulk_vafs_PL) <- bulk_keys[mask]

# Prepare result list

shared_stats <- lapply(names(branch_vcf), function(branch_name) {
  branch <- branch_vcf[[branch_name]]
  branch_keys <- variant_key(branch)
  
  # Get intersection
  shared <- intersect(branch_keys, names(bulk_vafs_PL))
  
  # Get VAFs of shared variants from bulk
  shared_vafs <- bulk_vafs_PL[shared]
  
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

group1 <- c("o") # BL trunk
group2 <- c("j") # BL intermediate PL
group3 <- c("c", "b", "e", "h", "g") # BL private PL
group4 <- c("p", "r") # BL intermediate BM
group5 <- c("J", "I", "L", "H", "O", "G", "R", "U", "T", "W", "Y", "C", "B", "A", "C2", "A2", "F", "u", "t", "w", "q") # BL private BM
group6 <- c("l", "k") # WT

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
  data.frame(VAF = vafs_group2, Group = "BL-Intermediate-PL"),
  data.frame(VAF = vafs_group3, Group = "BL-Private-PL"),
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
  cbind(get_shared_vs_total_stats(group2, shared_stats, tree_mut_counts), Group = "BL-Intermediate-PL"),
  cbind(get_shared_vs_total_stats(group3, shared_stats, tree_mut_counts), Group = "BL-Private-PL"),
  cbind(get_shared_vs_total_stats(group4, shared_stats, tree_mut_counts), Group = "BL-Intermediate-BM"),
  cbind(get_shared_vs_total_stats(group5, shared_stats, tree_mut_counts), Group = "BL-Private-BM"),
  cbind(get_shared_vs_total_stats(group6, shared_stats, tree_mut_counts), Group = "WT")
)

plot_df$Group <- factor(plot_df$Group, levels = c("BL-Trunk", "BL-Intermediate-PL", "BL-Private-PL", "BL-Intermediate-BM","BL-Private-BM", "WT"))
summary_stats$Group <- factor(summary_stats$Group, levels = c("BL-Trunk", "BL-Intermediate-PL", "BL-Private-PL", "BL-Intermediate-BM","BL-Private-BM", "WT"))

group_cols <- c(
  "BL-Trunk" = "#1D3557",
  "BL-Intermediate-PL" = "#0A9086",
  "BL-Private-PL" = "#A62639",
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
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_bulk_VAF_distribution_histogram_PL.pdf",
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
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_bulk_VAF_distribution_violin_LN.pdf",
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
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_bulk_VAF_distribution_intermediate_branches_violin_PL.pdf",
  plot = p6a,
  width = 8,
  height = 6
)

# now repeat for BM bulk

vaf_BM      <- as.numeric(vaf_mat[, 1])
mask      <- !is.na(vaf_BM) & vaf_BM >= 0.15

bulk_vafs_BM <- vaf_BM[mask]
names(bulk_vafs_BM) <- bulk_keys[mask]

# Prepare result list

shared_stats <- lapply(names(branch_vcf), function(branch_name) {
  branch <- branch_vcf[[branch_name]]
  branch_keys <- variant_key(branch)
  
  # Get intersection
  shared <- intersect(branch_keys, names(bulk_vafs_BM))
  
  # Get VAFs of shared variants from bulk
  shared_vafs <- bulk_vafs_BM[shared]
  
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
  #data.frame(VAF = vafs_group2, Group = "BL-Intermediate-PL") # no mutations
  data.frame(VAF = vafs_group3, Group = "BL-Private-PL"), 
  data.frame(VAF = vafs_group4, Group = "BL-Intermediate-BM"), 
  data.frame(VAF = vafs_group5, Group = "BL-Private-BM")
  #data.frame(VAF = vafs_group6, Group = "WT")  # no mutations
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
  cbind(get_shared_vs_total_stats(group2, shared_stats, tree_mut_counts), Group = "BL-Intermediate-PL"),
  cbind(get_shared_vs_total_stats(group3, shared_stats, tree_mut_counts), Group = "BL-Private-PL"),
  cbind(get_shared_vs_total_stats(group4, shared_stats, tree_mut_counts), Group = "BL-Intermediate-BM"),
  cbind(get_shared_vs_total_stats(group5, shared_stats, tree_mut_counts), Group = "BL-Private-BM"),
  cbind(get_shared_vs_total_stats(group6, shared_stats, tree_mut_counts), Group = "WT")
)

plot_df$Group <- factor(plot_df$Group, levels = c("BL-Trunk", "BL-Intermediate-PL", "BL-Private-PL", "BL-Intermediate-BM","BL-Private-BM", "WT"))
summary_stats$Group <- factor(summary_stats$Group, levels = c("BL-Trunk", "BL-Intermediate-PL", "BL-Private-PL", "BL-Intermediate-BM","BL-Private-BM", "WT"))

group_cols <- c(
  "BL-Trunk" = "#1D3557",
  "BL-Intermediate-PL" = "#0A9086",
  "BL-Private-PL" = "#A62639",
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
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_bulk_VAF_distribution_histogram_BM.pdf",
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
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_bulk_VAF_distribution_violin_BM.pdf",
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
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/P856/P856_bulk_VAF_distribution_intermediate_branches_violin_BM.pdf",
  plot = p6b,
  width = 8,
  height = 6
)

# Plot signatures on tree

# Get signatures 

all_signatures = get_known_signatures()
sbsblood <- read.table("~/Downloads/sigfit_cosmic3_bloodsig_Aug2020.txt", sep = "\t", header = T)
sbsblood = as.matrix(sbsblood)
SBSblood <- as.numeric(sbsblood[,"Signature.Blood"])
signatures <- cbind(SBSblood, all_signatures)

# Refitting part 

sub_sig <- signatures[, c("SBS1","SBS7a", "SBS9", "SBS17b", "SBS18", "SBSblood")] # from de novo extraction

contribution <- fit_to_signatures_strict_tree(mut_matrix = branch_mm, signatures = sub_sig, max_delta = 0.002, remove_min = 0)

pal <- c("#D2BD96","#6C5B78", "#0A9086", "#B3B3B3", "#A62639", "#1D3557")
names(pal) <- c("SBS1", "SBS7a","SBS9","SBS17b","SBS18","SBSblood")
sig_cols <- pal[rownames(contribution)]

plot_contribution(contribution, palette = pal)

# Check if all branches are explained well with these signatures: cosine < 0.85 with > 200 mutations would suggest you miss a mutation

p7 <- check_reconstructed_cosine_Alex(contribution, branch_mm, sub_sig, tree) +
  ylim(0,1)

# Save as PDF

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/P856/P856_reconstructed_cosine_per_branch.pdf",
  plot = p7,
  width = 8,
  height = 6
)

write.csv(contribution, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Trees/P856_per_branch_sbs_contribution.csv", row.names = TRUE)

# add signature contributions to your tree

tree = add_contribution(tree, contribution = contribution) # if you already did signature_fitting

p9 <- plot_tree_contribution_bars_new_alex(
  tree = tree,
  signatures = sub_sig,
  mut_matrix = branch_mm,
  signature_colors = sig_cols,
  title = "P856",
  bar_height = 0.02,
  scaling = 0.8,
  branch_color = "#000000",
  x_limits = c(0, 3600),
  x_breaks = seq(0, 3600, by = 200)
)

# Save as PDF

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/P856/P856_signatures_in_tree.pdf",
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
  column_names %in% group2 ~ "BL-Intermediate-PL",
  column_names %in% group3 ~ "BL-Private-PL",
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

# Re-order

contribution_wide <- contribution_wide |>
  dplyr::select(sample, WT, "BL-Trunk", "BL-Intermediate-PL", "BL-Private-PL", "BL-Intermediate-BM", "BL-Private-BM" )

# Convert the result into a matrix, using sample names as row names

contribution_matrix <- as.matrix(contribution_wide %>% column_to_rownames("sample"))

# View the resulting matrix

contribution_matrix

write.csv(contribution_matrix, "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/sbs_contr_per_group_P856.csv", row.names = T)

p10 <- plot_contribution(contribution_matrix[,], palette = sig_cols) 

# Save as PDF

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/P856/P856_signatures_per_group.pdf",
  plot = p10,
  width = 8,
  height = 6
)

# Convert to long format

df_long <- melt(contribution_matrix[,], varnames = c("Signature", "Sample"), value.name = "Contribution") %>%
  group_by(Sample) %>%
  mutate(Prop = Contribution / sum(Contribution)) %>%
  ungroup()

p11 <- ggplot(df_long, aes(x = Sample, y = Prop, fill = Signature)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = sig_cols) +
  scale_y_continuous(limits = c(0, 1)) +     # y-axis from 0 to 1
  labs(title = "Signature contribution per sample",
       x = "",
       y = "")

# Save as PDF

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Trees/P856/P856_signatures_per_group_separate_bars.pdf",
  plot = p11,
  width = 8,
  height = 6
)

