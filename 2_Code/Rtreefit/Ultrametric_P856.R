################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to generate ultrametric trees for Donor P856
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

tree = readRDS("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/P856/Filtered_samples_v5/CPW_01/TreeObject0.1.RDS")
#vcf = VariantAnnotation::readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/P856/Filtered_samples_v5/P856.vep.sub.SNV.autosomal.noBulksNoBlacklist.vcf")

# prepare tree

tree = prepare_tree(tree)

# plot bare tree

plot_gg_tree_base(tree)
plot_gg_tree_base(tree)+ geom_nodelab()
plot_gg_tree(tree, add_branch_length = TRUE, add_bootstrap = F, add_tip_label = F,add_title = "P856")

# VCF per branch
P856_contribution <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Trees/P856_per_branch_sbs_contribution.csv", row.names = 1)
P856_contribution$Signature <- rownames(P856_contribution)

subset_sum <- P856_contribution %>%
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
x[30,] #this row contains NAs (node 15-to-15 connection)
x<- x[-30,] #deleting NA-filled row
phylo_tree$edge.length<-unname(x$branch_length)
phylo_tree$edge
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label
) %>%
  mutate(age = case_when(
    grepl("LLP4", tip.label) ~ 14.9, #14.0 + 0.9
    grepl("DUB", tip.label) ~ 14.9, 
    TRUE ~ 15.6 # 14.7 +0.9
  ))
phylo_tree$agedf<-agedf
ggtree(phylo_tree)

#x <- x %>%
#  mutate(ctg_branch_lengths = ctg_totals[as.character(branch_id)])

x$sbs1_sbsblood_branch_lengths <- vals[as.character(x$branch_id)]


# correct for callable loci

# Step 1: Create sensitivity column in callable_df
max_callable <- 2745186691
input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') #dataframe
input_df_P856 <- input_df[input_df$Novogene_ID == "P856" & input_df$ResolveDNA_version %in% c("v1", "v2.0", "v2"),]
low_callable_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv')
below_curve_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')
filtered_samples <- unique(c(low_callable_df$Sample_name, below_curve_df$Sample_name, bad_baf_df$Sample_name , fail_vaf_df$samplename))  # samples that didn't pass QC
input_df_P856_filtered <- input_df_P856 %>% filter(!Sample_name %in% filtered_samples)
input_df_P856_filtered$Callable_Loci <- as.numeric(input_df_P856_filtered$Callable_Loci)
input_df_P856_filtered$sensitivity <- input_df_P856_filtered$Callable_Loci / max_callable

# Step 2: Define the product sensitivity calculator
calculate_product_sensitivity <- function(samples_str, sensitivity_df) {
  sample_list <- strsplit(samples_str, "\\|")[[1]]
  sensitivities <- sensitivity_df$sensitivity[sensitivity_df$Sample_name %in% sample_list]
  if (length(sensitivities) == 0) return(1)  # fallback to avoid division by zero
  return(1 - prod(1 - sensitivities, na.rm = TRUE))
}

# Step 3: Apply it to your tibble `x`
#x <- x %>%
#  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_P856_filtered),
#         corr_branch_lengths = ctg_branch_lengths / product_sensitivity)

x <- x %>%
  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_P856_filtered),
         corr_branch_lengths = sbs1_sbsblood_branch_lengths / product_sensitivity)

x <- x %>%
  mutate(corr_branch_lengths = round(corr_branch_lengths))

phylo_tree$edge.length<-unname(x$corr_branch_lengths)
phylo_tree$edge
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label
) %>%
  mutate(age = case_when(
    grepl("LLP4", tip.label) ~ 14.9, #14.0 + 0.9
    grepl("DUB", tip.label) ~ 14.9, 
    TRUE ~ 15.6 # 14.7 +0.9
  ))
phylo_tree$agedf<-agedf
phylo_tree$edge.length[is.na(phylo_tree$edge.length)] <- 0

p <- ggtree(phylo_tree) +
  geom_text2(aes(label = round(branch.length, 2)), hjust = -0.3) +
  theme_tree2()

print(p)

res_P856=fit_tree(tree=phylo_tree,switch_nodes = c(50,48),xcross = c(0.999,0.999), niter = 20000, cores = 4, model = "poisson_tree", early_growth_model_on = 1.0)

# Save ultrametric tree

# save(res_P856, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/res_P856.RData")
# load("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/res_P856.RData")

# Create the plot
p_P856 <- ggtree(res_P856$ultratree, color = "grey") +
  theme_tree2() +
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 16), breaks = 0:16) +
  # translucent grey box
  annotate(
    "rect",
    xmin = res_P856$upper_node_lims$lb95, xmax = res_P856$upper_node_lims$ub95,
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.2
  ) +
  # dotted vertical lines at x=0 and x=4.1
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 14.9, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 15.6, color = "red", linetype = "dotted") +
  geom_vline(xintercept = res_P856$upper_node_lims$mean, color = "black", linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
  )

print(p_P856)

res_P856$upper_node_lims$mean
res_P856$upper_node_lims$ub95
res_P856$upper_node_lims$lb95


# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/P856_ultrametric_tree_sbs1_blood_node50_48_0999.pdf",
  plot = p_P856,
  width = 8,
  height = 6
)

res_P856_v2=fit_tree(tree=phylo_tree,switch_nodes = c(48),xcross = c(0.999), niter = 20000, cores = 4, model = "poisson_tree", early_growth_model_on = 1.0)

# Save ultrametric tree

# save(res_P856, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/res_P856.RData")
# load("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/res_P856.RData")

# Create the plot
p_P856_v2 <- ggtree(res_P856_v2$ultratree, color = "grey") +
  theme_tree2() +
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 16), breaks = 0:16) +
  # translucent grey box
  annotate(
    "rect",
    xmin = res_P856_v2$upper_node_lims$lb95, xmax = res_P856_v2$upper_node_lims$ub95,
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.2
  ) +
  # dotted vertical lines at x=0 and x=4.1
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 14.9, color = "red", linetype = "dotted") +
  geom_vline(xintercept = 15.6, color = "red", linetype = "dotted") +
  geom_vline(xintercept = res_P856_v2$upper_node_lims$mean, color = "black", linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
  )

print(p_P856_v2)

res_P856_v2$upper_node_lims$mean
res_P856_v2$upper_node_lims$ub95
res_P856_v2$upper_node_lims$lb95


# # Save as PDF
ggsave(
filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/P856_ultrametric_tree_sbs1_blood_node_48_0999.pdf",
   plot = p_P856_v2,
   width = 8,
   height = 6
 )


# Calculate mutation rates for malignant and normal B-cells 

P856_malignant_cells_PL <- c( "PB14458-BLPL-BCELLP4L3", "PB14458-BLPL-BCELLP4K5", "PB14458-BLPL-BCELLP4M3", "PB14458-BLPL-BCELLP4B3", "PB14458-BLPL-BCELLP4L5" )

P856_malignant_cells_BM <- c( "P856GDDBBC63", "PB14458-BLBM-BCELLP2F2", "PB14458-BLBM-BCELLP2B3", "P856GDDBBC60", "P856GDDBBC46", "PB14458-BLBM-BCELLP2B4", "PB14458-BLBM-BCELLP2C4", "PB14458-BLBM-BCELLP2N2", "P856GDDBBC59", "P856GDDBBC64", "PB14458-BLBM-BCELLP2F4", "P856GDDBBC58", "P856GDDBBC54", "P856GDDBBC48", "PB14458-BLBM-BCELLP2I2", "PB14458-BLBM-BCELLP2L4",  "PB14458-BLBM-BCELLP2E4", "P856GDDBBC57", "P856GDDBBC62" ,"PB14458-BLBM-BCELLP2L3", "P856GDDBBC61")

P856_normal_cells <- c( "P856GDDUBC42", "P856GDDUBC44")

age_P856_diag <- 14.9
age_P856_rel <- 15.6

# Long-form expansion of `samples`, handling "," or "|" and trimming spaces
x_long <- x %>%
  mutate(row_id = row_number()) %>%
  separate_rows(samples, sep = "\\s*[|,]\\s*") %>%
  mutate(samples = trimws(samples)) %>%
  filter(samples != "")

# Keep only the malignant samples of LN
x_long_mal_PL <- x_long %>%
  filter(samples %in% P856_malignant_cells_PL)

# Sum corr_branch_lengths across ALL rows where each sample appears
per_sample_totals_PL <- x_long_mal_PL %>%
  group_by(samples) %>%
  summarise(
    total_corr_branch_length = sum(corr_branch_lengths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_corr_branch_length))

# MRCA length

MRCA_length_PL <- x %>%
  filter(branch_id %in% c("o", "j")) %>%
  summarise(MRCA = sum(corr_branch_lengths, na.rm = TRUE)) %>%
  pull(MRCA)

# Get mutations after MRCA

per_sample_totals_PL$MRCA <- MRCA_length_PL

# Remove MRCA mutation load per sample
per_sample_totals_PL <- per_sample_totals_PL %>%
  mutate(total_corr_branch_length_wo_MRCA = total_corr_branch_length - MRCA)

# Compute latency time
latency_time_P856_PL <- age_P856_diag - 6.100166

# Compute mutation rate *per malignant sample*
per_sample_totals_PL <- per_sample_totals_PL %>%
  mutate(mutation_rate_post_expansion = total_corr_branch_length_wo_MRCA / latency_time_P856_PL)

# Get mean and SE
mean_mutation_rate_post_expansion_P856_PL <- mean(per_sample_totals_PL$mutation_rate_post_expansion, na.rm = TRUE)
se_mutation_rate_post_expansion_P856_PL <- sd(per_sample_totals_PL$mutation_rate_post_expansion, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals_PL$mutation_rate_post_expansion)))

mean_mutation_rate_post_expansion_P856_PL
se_mutation_rate_post_expansion_P856_PL

# Mutation rate of pre-expansion cell
mutation_rate_pre_expansion_P856_PL <- MRCA_length_PL/6.100166
mutation_rate_pre_expansion_P856_PL

# Keep only the malignant samples of BM
x_long_mal_BM <- x_long %>%
  filter(samples %in% P856_malignant_cells_BM)

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
  filter(branch_id %in% c("o", "p")) %>%
  summarise(MRCA = sum(corr_branch_lengths, na.rm = TRUE)) %>%
  pull(MRCA)

# Get mutations after MRCA

per_sample_totals_BM$MRCA <- MRCA_length_BM

# Remove MRCA mutation load per sample
per_sample_totals_BM <- per_sample_totals_BM %>%
  mutate(total_corr_branch_length_wo_MRCA = total_corr_branch_length - MRCA)

# Compute latency time
latency_time_P856_BM <- age_P856_rel - res_P856_v2$upper_node_lims$mean

# Compute mutation rate *per malignant sample*
per_sample_totals_BM <- per_sample_totals_BM %>%
  mutate(mutation_rate_post_expansion = total_corr_branch_length_wo_MRCA / latency_time_P856_BM)

# Get mean and SE
mean_mutation_rate_post_expansion_P856_BM <- mean(per_sample_totals_BM$mutation_rate_post_expansion, na.rm = TRUE)
se_mutation_rate_post_expansion_P856_BM <- sd(per_sample_totals_BM$mutation_rate_post_expansion, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals_BM$mutation_rate_post_expansion)))

mean_mutation_rate_post_expansion_P856_BM
se_mutation_rate_post_expansion_P856_BM

# Mutation rate of pre-expansion cell
mutation_rate_pre_expansion_P856_BM <- MRCA_length_BM/res_P856_v2$upper_node_lims$mean
mutation_rate_pre_expansion_P856_BM

# Keep only the normal samples 
x_long_norm <- x_long %>%
  filter(samples %in% P856_normal_cells)

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
  mutate(mutation_rate = total_corr_branch_length / age_P856_diag)

# Get mean and SE
mean_mutation_rate_normal_P856 <- mean(per_sample_totals_norm$mutation_rate, na.rm = TRUE)
se_mutation_rate_normal_P856 <- sd(per_sample_totals_norm$mutation_rate, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals_norm$mutation_rate)))

mean_mutation_rate_normal_P856
se_mutation_rate_normal_P856


mean(per_sample_totals_BM$total_corr_branch_length_wo_MRCA)
mean(per_sample_totals_PL$total_corr_branch_length_wo_MRCA)
mean(per_sample_totals_BM$total_corr_branch_length_wo_MRCA) - mean(per_sample_totals_PL$total_corr_branch_length_wo_MRCA)
(mean(per_sample_totals_BM$total_corr_branch_length_wo_MRCA) - mean(per_sample_totals_PL$total_corr_branch_length_wo_MRCA))/0.7
