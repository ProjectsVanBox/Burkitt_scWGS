################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at CellPhy Tree of Donor P3G6 (branch VCFs for manual checking, choosing percentage, signatures per branch)
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
library(stringdist)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(ggnewscale)
ref_genome = 'BSgenome.Hsapiens.NCBI.GRCh38'
library(MutationalPatterns)
library(VariantAnnotation)
library(ggplot2)
library(dplyr)
library(cellPhyWrapperPlotting) #devtools::install_local("~/surfdrive/Shared/pmc_vanboxtel/general/2_Bioinformatics/Scripts/cellPhyWrapperPlotting/",force = TRUE)
source("~/surfdrive/Shared/pmc_vanboxtel/personal/asteemers/R_packages/plot_signature_contribution_new.R")
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/PTATO/GeneralFunctions.R')
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')
theme_set(theme_classic())
source("~/surfdrive/Shared/pmc_vanboxtel/general/2_Bioinformatics/colors/Jurrians_colors.R")

# Load tree object, vcf used for tree building and bulk vcf

tree = readRDS("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/P3G6/Filtered_samples/CPW_04/TreeObject0.4.RDS")
vcf = VariantAnnotation::readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/P3G6/P3G6.vep.sub.SNV.autosomal.noBulksNoBlacklist.vcf")
vcf_bulk = VariantAnnotation::readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P3G6/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P3G6_bulk.vep.SMuRF.filtered.sorted.vcf.gz")

# prepare tree
tree = prepare_tree(tree)

# plot bare tree
plot_gg_tree_base(tree)
plot_gg_tree(tree, add_branch_length = TRUE, add_bootstrap = F, add_tip_label = F,add_title = "PRN4")

# VCF per branch
branch_vcf = extract_vcf_per_branch(tree = tree, vcf = vcf, ref_genome = ref_genome)
branch_grl = convert_vcf_to_granges(branch_vcf_list = branch_vcf, ref_genome = ref_genome)
branch_mm = mut_matrix(branch_grl, ref_genome)

# Step 1: save branch vcfs to manually check "odd" cells i.e. ones that branch off earlier than the major clone

# In this case I wan to check branch h for PRN4GPDLBC17 and PRN4GPDLBC15 and branch e for PB08410-BLBM-BCELLP2G8
# Specifically check for linked SNPs and badly covered regions

lapply(c("h", "e"), function(name) {
  filename <- file.path("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PRN4/Filtered_samples_v2/branch_vcfs_to_manually_check", paste0(name, "_branch.vcf"))
  VariantAnnotation::writeVcf(branch_vcf[[name]], file = filename)
})


# Preprocess bulk samples

vcf_bulk_snvs <- vcf_bulk[isSNV(vcf_bulk)] # only SNVs
autosomes <- as.character(1:22) 
vcf_bulk_autosomal <- vcf_bulk_snvs[seqnames(rowRanges(vcf_bulk_snvs)) %in% autosomes] # only autosomes 
vaf_matrix <- geno(vcf_bulk_autosomal)$VAF
vaf_values <- as.numeric(vaf_matrix[, 1])
vcf_bulk_filtered <- vcf_bulk_autosomal[vaf_values > 0.3] # filter on VAF (cutoff decided in previous script now called Mutational_load_VAF_cutoff)

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

bulk_keys <- variant_key(vcf_bulk_filtered)
bulk_vaf_matrix <- geno(vcf_bulk_filtered)$VAF
bulk_vafs <- as.numeric(bulk_vaf_matrix[, 1]) 

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

group1 <- c("E3")
group2 <- c("F3", "I3","J2","N","Y", "Z","W","I2", "H2","G2")
group3 <- c("A", "G3", "H3", "X", "O", "P", "U", "R", "S", "K","L","B","H","F","C","D","a","b","c","d","i","f","g","D2","A2","B2","I","m","u","v","s","q","n","o")

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

# Build one long data frame for ggplot
plot_df <- rbind(
  data.frame(VAF = vafs_group1, Group = "Trunk"),
  data.frame(VAF = vafs_group2, Group = "Intermediate"),
  data.frame(VAF = vafs_group3, Group = "Private")
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
  cbind(get_shared_vs_total_stats(group1, shared_stats, tree_mut_counts), Group = "Trunk"),
  cbind(get_shared_vs_total_stats(group2, shared_stats, tree_mut_counts), Group = "Intermediate"),
  cbind(get_shared_vs_total_stats(group3, shared_stats, tree_mut_counts), Group = "Private")
)

plot_df$Group <- factor(plot_df$Group, levels = c("Trunk", "Intermediate", "Private"))
summary_stats$Group <- factor(summary_stats$Group, levels = c("Trunk", "Intermediate", "Private"))


ggplot(plot_df, aes(x = VAF)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  facet_wrap(~Group, scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  geom_text(
    data = summary_stats,
    aes(x = 0.8, y = Inf, label = Label),
    inherit.aes = FALSE,
    vjust = 2, size = 4.5, fontface = "bold"
  ) +
  theme_classic() +
  labs(
    title = "Bulk VAF Distributions of Shared Mutations",
    x = "Bulk VAF",
    y = "Count"
  )
tree = readRDS("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Markus/P3G6/CPW_01/TreeObject0.1.RDS")
vcf = VariantAnnotation::readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Markus/P3G6/ptato_asap_woBCELLBULK.vcf.gz")

# prepare tree
tree = prepare_tree(tree)

# plot bare tree
plot_gg_tree_base(tree)
plot_gg_tree(tree, add_branch_length = TRUE, add_bootstrap = F, add_tip_label = F,add_title = "P3G6 | Ascites | 13.7Y")

# VCF per branch
branch_vcf = extract_vcf_per_branch(tree = tree, vcf = vcf, ref_genome = ref_genome)
branch_grl = convert_vcf_to_granges(branch_vcf_list = branch_vcf, ref_genome = ref_genome)

# Function to extract VAFs from a single VCF
extract_vaf <- function(vcf_per_branch) {
  AD_values <- vcf_per_branch@assays@data$AD  # Extract Allelic Depth
  
  # Convert AD to numeric matrix if it's a list
  if (is.list(AD_values)) {
    AD_values <- lapply(AD_values, function(x) {
      if (is.null(x) || length(x) < 2 || !all(sapply(x, is.numeric))) {
        return(NA)  # Handle missing or malformed values
      }
      return(as.numeric(x))  # Ensure numeric conversion
    })
    AD_values <- do.call(cbind, AD_values)  # Convert list to matrix
  }
  
  # Ensure it's a numeric matrix and has at least two rows
  if (!is.matrix(AD_values) || nrow(AD_values) < 2) {
    return(NA)  # Return NA if the structure is invalid
  }
  
  # Calculate VAF safely
  vaf_values <- apply(AD_values, 2, function(x) {
    if (length(x) < 2 || sum(x, na.rm = TRUE) == 0) {
      return(NA)  # Avoid division by zero
    }
    return(as.numeric(x[2]) / sum(as.numeric(x), na.rm = TRUE))
  })
  
  vaf_values <- vaf_values[!is.na(vaf_values) & vaf_values != 0.0]
}

# profile per branch --> can be used for any MutationalPatterns analysis separate from the tree
branch_mm = mut_matrix(branch_grl, ref_genome)


### Get signatures 
signatures = get_known_signatures()
#pta_v1_sig = read.table("~/surfdrive - Alexander Steemers@surfdrive.surf.nl/Shared/pmc_vanboxtel/projects/Burkitt_lymphoma/1_Input/WGS/PTA_Artefact_Signature.txt", sep = "\t", header = T)
#pta_v1_sig = as.matrix(pta_v1_sig)
#PTA_v1 <- as.numeric(pta_v1_sig[,"PTA"])
#PTA_v1 <- PTA_v1[!is.na(PTA_v1)]
#pta_v2_sig = read.table("~/hpc/pmc_vanboxtel/resources/signatures/PTAv2_Artefact_Signature.txt", sep = "\t", header = T) 
#pta_v2_sig = as.matrix(pta_v2_sig) 
#PTA_v2 <- as.numeric(pta_v2_sig[,"PTAv2"])
#hspc_sig = read.table("~/surfdrive - Alexander Steemers@surfdrive.surf.nl/Shared/pmc_vanboxtel/projects/Burkitt_lymphoma/1_Input/WGS/sigProfiler_SBS_working_signatures_incl_hspc.txt", sep = "\t", header = T)
#hspc_sig = as.matrix(hspc_sig)
#HSPC <- as.numeric(hspc_sig[,"HSPC"])
sbsblood <- read.table("~/Downloads/sigfit_cosmic3_bloodsig_Aug2020.txt", sep = "\t", header = T)
sbsblood = as.matrix(sbsblood)
SBSblood <- as.numeric(sbsblood[,"Signature.Blood"])
signatures <- cbind(SBSblood, signatures)

### Refitting part 
sub_sig <- signatures[, c("SBSblood", "SBS1", "SBS7a","SBS8", "SBS9", "SBS17b", "SBS18")] #signatures from Machado paper
sub_sig <- signatures[, c("SBSblood", "SBS1","SBS8", "SBS9")] #signatures from Machado paper

contribution <- fit_to_signatures_strict_tree(mut_matrix = branch_mm, signatures = sub_sig, max_delta = 0.01, remove_min = 20)
contribution <- fit_to_signatures_strict_tree(mut_matrix = branch_mm, signatures = sub_sig, max_delta = 0.01, remove_min = 0)
write.csv(contribution, "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_lymphoma/3_Output/CellPhy_mutationalpatterns/sbs_contr_per_branch_P3G6.csv", row.names = T)
plot_contribution(contribution, palette = dist_cols50)

# check if all branches are explained well with these signatures: cosine < 0.85 with > 200 mutations would suggest you miss a mutation
check_reconstructed_cosine(contribution, branch_mm, sub_sig, tree) +
  ylim(0,1)

# add signature contributions to your tree
tree = add_contribution(tree, contribution = contribution) # if you already did signature_fitting

# plot bars with all signatures in the branches of the tree
plot_tree_contribution_bars_new(tree = tree, signatures = sub_sig, mut_matrix = branch_mm, title = "P3G6 | Ascites | 13.7Y") 