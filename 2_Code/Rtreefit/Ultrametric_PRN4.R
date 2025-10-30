################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to generate ultrametric trees for Donor PRN4
# Author: Alexander Steemers
################################################################################

# Load libraries and functions

library(tidyverse)
library(vcfR) 
library(readxl) 
library(ggtree) 
library(treeio) 
library(ape)
library(rtreefit)
library(stringi)
library(stringdist) 
library(dplyr)
library(BSgenome.Hsapiens.NCBI.GRCh38) 
library(ggnewscale) 
ref_genome = 'BSgenome.Hsapiens.NCBI.GRCh38' 
library(MutationalPatterns) 
library(VariantAnnotation)
theme_set(theme_classic())
library(cellPhyWrapperPlotting)

# Load tree object

tree = readRDS("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PRN4/Filtered_samples/CPW_01/TreeObject0.1.RDS")
#vcf = VariantAnnotation::readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PRN4/PRN4.vep.sub.SNV.autosomal.noBulksNoBlacklist.vcf")

# prepare tree

tree = prepare_tree(tree)

# plot bare tree

plot_gg_tree_base(tree)
plot_gg_tree_base(tree)+ geom_nodelab()
plot_gg_tree(tree, add_branch_length = TRUE, add_bootstrap = F, add_tip_label = F,add_title = "PRN4")

# VCF per branch
PRN4_contribution <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Trees/PRN4_per_branch_sbs_contribution.csv", row.names = 1)
PRN4_contribution$Signature <- rownames(PRN4_contribution)

subset_sum <- PRN4_contribution %>%
  filter(Signature %in% c("SBS1", "SBSblood")) %>%
  summarise(across(where(is.numeric), sum))

vals <- setNames(as.numeric(subset_sum[1, ]), colnames(subset_sum))

#branch_vcf = extract_vcf_per_branch(tree = tree, vcf = vcf, ref_genome = ref_genome)
#branch_grl = convert_vcf_to_granges(branch_vcf_list = branch_vcf, ref_genome = ref_genome)
#branch_mm = mut_matrix(branch_grl, ref_genome)
#branch_mm_ctg <- branch_mm[grep("\\[C>T\\]G", rownames(branch_mm)), ]
#ctg_totals <- colSums(branch_mm_ctg)


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
  age = rep(5.0 ,52)  # Replace these with your actual ages (4.1 yo + 0.9)
)
phylo_tree$agedf<-agedf
ggtree(phylo_tree)


#x <- x %>%
#  mutate(ctg_branch_lengths = ctg_totals[as.character(branch_id)])

x$sbs1_sbsblood_branch_lengths <- vals[as.character(x$branch_id)]


# correct for callable loci

# Step 1: Create sensitivity column in callable_df
max_callable <- 2745186691
input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') #dataframe
input_df_PRN4 <- input_df[input_df$Novogene_ID == "PRN4" & input_df$ResolveDNA_version %in% c("v1", "v2.0", "v2"),]
low_callable_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv')
below_curve_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')
filtered_samples <- unique(c(low_callable_df$Sample_name, below_curve_df$Sample_name, bad_baf_df$Sample_name , fail_vaf_df$samplename))  # samples that didn't pass QC
input_df_PRN4_filtered <- input_df_PRN4 %>% filter(!Sample_name %in% filtered_samples)
input_df_PRN4_filtered$Callable_Loci <- as.numeric(input_df_PRN4_filtered$Callable_Loci)
input_df_PRN4_filtered$sensitivity <- input_df_PRN4_filtered$Callable_Loci / max_callable

# Step 2: Define the product sensitivity calculator
calculate_product_sensitivity <- function(samples_str, sensitivity_df) {
  sample_list <- strsplit(samples_str, "\\|")[[1]]
  sensitivities <- sensitivity_df$sensitivity[sensitivity_df$Sample_name %in% sample_list]
  if (length(sensitivities) == 0) return(1)  # fallback to avoid division by zero
  return(1 - prod(1 - sensitivities, na.rm = TRUE))
}

# Step 3: Apply it to your tibble `x`
#x <- x %>%
#  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_PRN4_filtered),
#         corr_branch_lengths = ctg_branch_lengths / product_sensitivity)

x <- x %>%
  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_PRN4_filtered),
         corr_branch_lengths = sbs1_sbsblood_branch_lengths / product_sensitivity)

x <- x %>%
  mutate(corr_branch_lengths = round(corr_branch_lengths))

phylo_tree$edge.length<-unname(x$corr_branch_lengths)
phylo_tree$edge
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label,  # Must match tree$tip.label
  age = rep(5.0 ,52)  # Replace these with your actual ages (4.1 yo + 0.9)
)
phylo_tree$agedf<-agedf
phylo_tree$edge.length[is.na(phylo_tree$edge.length)] <- 0

p <- ggtree(phylo_tree) +
  geom_text2(aes(label = round(branch.length, 2)), hjust = -0.3) +
  theme_tree2()

print(p)

res_PRN4=fit_tree(tree=phylo_tree,switch_nodes = c(43,48),xcross = c(0.999,0.999), niter = 20000, cores = 4, model = "poisson_tree", early_growth_model_on = 1.0)

# Save ultrametric tree

#save(res_PRN4, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/res_PRN4.RData")
#load("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/res_PRN4.RData")

# Create the plot
p_PRN4 <- ggtree(res_PRN4$ultratree, color = "grey") +
  theme_tree2() +
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 5.1), breaks = 0:5) +
  # translucent grey box
  annotate(
    "rect",
    xmin = res_PRN4$upper_node_lims$lb95, xmax = res_PRN4$upper_node_lims$ub95,
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.2
  ) +
  # dotted vertical lines at x=0 and x=4.1
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 5.0, color = "red", linetype = "dotted") +
  geom_vline(xintercept = res_PRN4$upper_node_lims$mean, color = "black", linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
  )

print(p_PRN4)


res_PRN4$upper_node_lims$mean
res_PRN4$upper_node_lims$ub95
res_PRN4$upper_node_lims$lb95

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/PRN4_ultrametric_tree_sbs1_blood_node43_48_0999.pdf",
  plot = p_PRN4,
  width = 8,
  height = 6
)


# Test: does removing a sample change the MRCA estimate?

tips_to_remove <- c("PRN4GPDLBC22")

phylo_tree_subset_22 <- drop.tip(phylo_tree, tips_to_remove )

p <- ggtree(phylo_tree_subset_22) +
  geom_text2(aes(label = round(branch.length, 2)), hjust = -0.3) +
  theme_tree2()

print(p)

res_PRN4_22=fit_tree(tree=phylo_tree_subset_22,switch_nodes = c(31),xcross = c(0.999), niter = 1000, cores = 4, model = "poisson_tree", early_growth_model_on = 1.0)

# Create the plot
p22 <- ggtree(res_PRN4_22$ultratree, color = "grey") +
  theme_tree2() +
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 5.1), breaks = 0:5) +
  # translucent grey box
  annotate(
    "rect",
    xmin = res_PRN4_22$upper_node_lims$lb95, xmax = res_PRN4_22$upper_node_lims$ub95,
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.2
  ) +
  # dotted vertical lines at x=0 and x=4.1
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 5.0, color = "red", linetype = "dotted") +
  geom_vline(xintercept = res_PRN4_22$upper_node_lims$mean, color = "black", linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
  )

print(p22)

tips_to_remove <- c("PRN4GPDLBC19")

phylo_tree_subset_19 <- drop.tip(phylo_tree, tips_to_remove )

p <- ggtree(phylo_tree_subset_19) +
  geom_text2(aes(label = round(branch.length, 2)), hjust = -0.3) +
  theme_tree2()

print(p)

res_PRN4_19=fit_tree(tree=phylo_tree_subset_19,switch_nodes = c(31),xcross = c(0.999), niter = 1000, cores = 4, model = "poisson_tree", early_growth_model_on = 1.0)

# Create the plot
p19 <- ggtree(res_PRN4_19$ultratree, color = "grey") +
  theme_tree2() +
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 5.1), breaks = 0:5) +
  # translucent grey box
  annotate(
    "rect",
    xmin = res_PRN4_19$upper_node_lims$lb95, xmax = res_PRN4_19$upper_node_lims$ub95,
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.2
  ) +
  # dotted vertical lines at x=0 and x=4.1
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 5.0, color = "red", linetype = "dotted") +
  geom_vline(xintercept = res_PRN4_19$upper_node_lims$mean, color = "black", linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
  )

print(p19)


res_PRN4$upper_node_lims$mean
res_PRN4_22$upper_node_lims$mean
res_PRN4_19$upper_node_lims$meam

# Save PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/PRN4_ultrametric_tree_sbs1_blood_node45_48_43_0999.pdf",
  plot = p,
  width = 8,
  height = 6
)

# Calculate mutation rates for malignant and normal B-cells 

PNR4_malignant_cells_LN <- c( "PRN4GPDLBC15", "PRN4GPDLBC17","PRN4GPDLBC09","PRN4GPDLBC11", "PB08410-BLLN-BCELLP4F10",  "PRN4GPDLBC21", "PRN4GPDLBC23", "PB08410-BLLN-BCELLP4D10", "PB08410-BLLN-BCELLP2D10","PRN4GPDLBC22", "PB08410-BLLN-BCELLP4B11", "PB08410-BLLN-BCELLP1B11", "PRN4GPDLBC20", "PB08410-BLLN-BCELLP2B10","PB08410-BLLN-BCELLP2E10",  "PRN4GPDLBC10", "PB08410-BLLN-BCELLP4G10", "PRN4GPDLBC16", "PRN4GPDLBC19")

PNR4_malignant_cells_BM <- c("PB08410-BLBM-BCELLP2G8", "PB08410-BLBM-BCELLP2D8", "PB08410-BLBM-BCELLP5E8", "PB08410-BLBM-BCELLP5F8")
PNR4_normal_cells <- c( "PB08410-BLBM-BCELLP5G10", "PRN4GPDBBC07", "PB08410-BLBM-BCELLP2C8")

age_PNR4 <- 5.0

# Long-form expansion of `samples`, handling "," or "|" and trimming spaces
x_long <- x %>%
  mutate(row_id = row_number()) %>%
  separate_rows(samples, sep = "\\s*[|,]\\s*") %>%
  mutate(samples = trimws(samples)) %>%
  filter(samples != "")

# Keep only the malignant samples of LN
x_long_mal_LN <- x_long %>%
  filter(samples %in% PNR4_malignant_cells_LN)

# Sum corr_branch_lengths across ALL rows where each sample appears
per_sample_totals_LN <- x_long_mal_LN %>%
  group_by(samples) %>%
  summarise(
    total_corr_branch_length = sum(corr_branch_lengths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_corr_branch_length))

# MRCA length

MRCA_length_LN <- x %>%
  filter(branch_id %in% c("h", "p")) %>%
  summarise(MRCA = sum(corr_branch_lengths, na.rm = TRUE)) %>%
  pull(MRCA)

# Get mutations after MRCA

per_sample_totals_LN$MRCA <- MRCA_length_LN

# Remove MRCA mutation load per sample
per_sample_totals_LN <- per_sample_totals_LN %>%
  mutate(total_corr_branch_length_wo_MRCA = total_corr_branch_length - MRCA)

# Compute latency time
latency_time_PNR4_LN <- age_PNR4 - res_PRN4$upper_node_lims$mean[1]

# Compute mutation rate *per malignant sample*
per_sample_totals_LN <- per_sample_totals_LN %>%
  mutate(mutation_rate_post_expansion = total_corr_branch_length_wo_MRCA / latency_time_PNR4_LN)

# Get mean and SE
mean_mutation_rate_post_expansion_PNR4_LN <- mean(per_sample_totals_LN$mutation_rate_post_expansion, na.rm = TRUE)
se_mutation_rate_post_expansion_PNR4_LN <- sd(per_sample_totals_LN$mutation_rate_post_expansion, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals_LN$mutation_rate_post_expansion)))

mean_mutation_rate_post_expansion_PNR4_LN
se_mutation_rate_post_expansion_PNR4_LN

# Mutation rate of pre-expansion cell
mutation_rate_pre_expansion_PNR4_LN <- MRCA_length_LN/res_PRN4$upper_node_lims$mean[1]
mutation_rate_pre_expansion_PNR4_LN

# Keep only the malignant samples of LN
x_long_mal_BM <- x_long %>%
  filter(samples %in% PNR4_malignant_cells_BM)

# Sum corr_branch_lengths across ALL rows where each sample appears
per_sample_totals_BM <- x_long_mal_BM %>%
  group_by(samples) %>%
  summarise(
    total_corr_branch_length = sum(corr_branch_lengths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_corr_branch_length))

# MRCA length

MRCA_length_BM <- x %>%
  filter(branch_id %in% c("h", "o")) %>%
  summarise(MRCA = sum(corr_branch_lengths, na.rm = TRUE)) %>%
  pull(MRCA)

# Get mutations after MRCA

per_sample_totals_BM$MRCA <- MRCA_length_BM

# Remove MRCA mutation load per sample
per_sample_totals_BM <- per_sample_totals_BM %>%
  mutate(total_corr_branch_length_wo_MRCA = total_corr_branch_length - MRCA)

# Compute latency time
latency_time_PNR4_BM <- age_PNR4 - res_PRN4$upper_node_lims$mean[2]

# Compute mutation rate *per malignant sample*
per_sample_totals_BM <- per_sample_totals_BM %>%
  mutate(mutation_rate_post_expansion = total_corr_branch_length_wo_MRCA / latency_time_PNR4_BM)

# Get mean and SE
mean_mutation_rate_post_expansion_PNR4_BM <- mean(per_sample_totals_BM$mutation_rate_post_expansion, na.rm = TRUE)
se_mutation_rate_post_expansion_PNR4_BM <- sd(per_sample_totals_BM$mutation_rate_post_expansion, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals_BM$mutation_rate_post_expansion)))

mean_mutation_rate_post_expansion_PNR4_BM
se_mutation_rate_post_expansion_PNR4_BM

# Mutation rate of pre-expansion cell
mutation_rate_pre_expansion_PNR4_BM <- MRCA_length_BM/res_PRN4$upper_node_lims$mean[2]
mutation_rate_pre_expansion_PNR4_BM

# Keep only the normal samples 
x_long_norm <- x_long %>%
  filter(samples %in% PNR4_normal_cells)

# Sum corr_branch_lengths across ALL rows where each sample appears
per_sample_totals_norm <- x_long_norm %>%
  group_by(samples) %>%
  summarise(
    total_corr_branch_length = sum(corr_branch_lengths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_corr_branch_length))

# Calculate mutation rate per normal sample
per_sample_totals_norm <- per_sample_totals_norm %>%
  mutate(mutation_rate = total_corr_branch_length / age_PNR4)

# Get mean and SE
mean_mutation_rate_normal_PNR4 <- mean(per_sample_totals_norm$mutation_rate, na.rm = TRUE)
se_mutation_rate_normal_PNR4 <- sd(per_sample_totals_norm$mutation_rate, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals_norm$mutation_rate)))

mean_mutation_rate_normal_PNR4
se_mutation_rate_normal_PNR4

