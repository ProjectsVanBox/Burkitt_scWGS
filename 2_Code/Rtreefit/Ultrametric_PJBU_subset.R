################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to generate ultrametric trees for Donor PJBU
# Author: Alexander Steemers
################################################################################

# Load libraries and functions
#.libPaths("/hpc/local/Rocky8/pmc_vanboxtel/bin/R-4.4.3/lib64/R/library")

library(tidyverse)
#library(vcfR) 
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

tree = readRDS("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PJBU/Filtered_samples_v3/CPW_01/TreeObject0.1.RDS")
#vcf = VariantAnnotation::readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PJBU/Filtered_samples_v3/PJBU.vep.sub.SNV.autosomal.noBulksNoBlacklist.vcf")

# prepare tree

tree = prepare_tree(tree)

# plot bare tree

plot_gg_tree_base(tree)
plot_gg_tree_base(tree)+ geom_nodelab()
plot_gg_tree(tree, add_branch_length = TRUE, add_bootstrap = F, add_tip_label = F,add_title = "PJBU")

# VCF per branch
PJBU_contribution <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Trees/PJBU_per_branch_sbs_contribution.csv", row.names = 1)
PJBU_contribution$Signature <- rownames(PJBU_contribution)

subset_sum <- PJBU_contribution %>%
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
x[60,] #this row contains NAs (node 15-to-15 connection)
x<- x[-60,] #deleting NA-filled row
phylo_tree$edge.length<-unname(x$branch_length)
phylo_tree$edge
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label,  # Must match tree$tip.label
  age = rep(7.5,116)  # Replace these with your actual ages 6.6 +0.9
)
phylo_tree$agedf<-agedf
ggtree(phylo_tree)

#x <- x %>%
#  mutate(ctg_branch_lengths = ctg_totals[as.character(branch_id)])

x$sbs1_sbsblood_branch_lengths <- vals[as.character(x$branch_id)]


# correct for callable loci

# Step 1: Create sensitivity column in callable_df
max_callable <- 2745186691
input_df <-  read_excel('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/1_Input/Sample_overview.xlsx') #dataframe
input_df_PJBU <- input_df[input_df$Novogene_ID == "PJBU" & input_df$ResolveDNA_version %in% c("v1", "v2.0", "v2"),]
low_callable_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/low_callable_loci.csv')
below_curve_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')
filtered_samples <- unique(c(low_callable_df$Sample_name, below_curve_df$Sample_name, bad_baf_df$Sample_name , fail_vaf_df$samplename))  # samples that didn't pass QC
input_df_PJBU_filtered <- input_df_PJBU %>% filter(!Sample_name %in% filtered_samples)
input_df_PJBU_filtered$Callable_Loci <- as.numeric(input_df_PJBU_filtered$Callable_Loci)
input_df_PJBU_filtered$sensitivity <- input_df_PJBU_filtered$Callable_Loci / max_callable

# Step 2: Define the product sensitivity calculator
calculate_product_sensitivity <- function(samples_str, sensitivity_df) {
  sample_list <- strsplit(samples_str, "\\|")[[1]]
  sensitivities <- sensitivity_df$sensitivity[sensitivity_df$Sample_name %in% sample_list]
  if (length(sensitivities) == 0) return(1)  # fallback to avoid division by zero
  return(1 - prod(1 - sensitivities, na.rm = TRUE))
}

# Step 3: Apply it to your tibble `x`
#x <- x %>%
#  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_PJBU_filtered),
#         corr_branch_lengths = ctg_branch_lengths / product_sensitivity)

x <- x %>%
  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_PJBU_filtered),
         corr_branch_lengths = sbs1_sbsblood_branch_lengths / product_sensitivity)

x <- x %>%
  mutate(corr_branch_lengths = round(corr_branch_lengths))

phylo_tree$edge.length<-unname(x$corr_branch_lengths)
phylo_tree$edge
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label,  # Must match tree$tip.label
  age = rep(7.5,116)  # Replace these with your actual ages 6.6 +0.9
)
phylo_tree$agedf<-agedf
phylo_tree$edge.length[is.na(phylo_tree$edge.length)] <- 0

plot(phylo_tree)
p <- ggtree(phylo_tree) +
  geom_nodelab() +
  geom_text2(aes(label = round(branch.length, 2)), hjust = -0.3) +
  theme_tree2()

print(p)

p <- ggtree(phylo_tree) +
  #geom_nodelab() +
  geom_text2(aes(label = round(branch.length, 2)), hjust = -0.3) +
  theme_tree2()

print(p)

tips_to_remove <- c("PJBUGTDBBC73") #remove one normal B-cell

phylo_tree_subset_73 <- drop.tip(phylo_tree, tips_to_remove )



# nodes_to_remove <- as.vector(unlist(x[which(x$tip.label %in% tips_to_remove),2]))
# nodes_to_remove <- c(nodes_to_remove,as.vector(unlist(x[which(x$branch_id %in% c("l","o")),"node"])))
# x2 <- x[-which(x$tip.label %in% tips_to_remove),]
# x2 <- x2[-which(x2$branch_id %in% c("l","o")),]
# 
# phylo_tree_subset$edge <- x2[,c(1:2)]
# phylo_tree_subset$edge <- phylo_tree_subset$edge[-which(phylo_tree_subset$edge$parent %in% nodes_to_remove),]

res_PJBU_73=fit_tree(tree=phylo_tree_subset_73,switch_nodes = c(67),xcross = c(0.999), niter = 20000, cores = 12, model = "poisson_tree", early_growth_model_on = 1.0)

# Create the plot
p73 <- ggtree(res_PJBU_73$ultratree, color = "grey") +
  theme_tree2() +
  #geom_tiplab()+
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 9), breaks = 0:9) +
  # translucent grey box
  annotate(
    "rect",
    xmin = res_PJBU_73$upper_node_lims$lb95, xmax = res_PJBU_73$upper_node_lims$ub95,
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.2
  ) +
  # dotted vertical lines at x=0 and x=4.1
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 7.5, color = "red", linetype = "dotted") +
  geom_vline(xintercept = res_PJBU_73$upper_node_lims$mean, color = "black", linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
  )

print(p73)

res_PJBU_73$upper_node_lims$mean
res_PJBU_73$upper_node_lims$ub95
res_PJBU_73$upper_node_lims$lb95

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/PJBU_ultrametric_tree_sbs1_blood_node67_0999_wo_PJBUGTDABC73.pdf",
  plot = p73,
  width = 8,
  height = 6
)

# tips_to_remove <- c("PJBUGTDABC25")
# 
# phylo_tree_subset_25 <- drop.tip(phylo_tree, tips_to_remove )
# 
# res_PJBU_25=fit_tree(tree=phylo_tree_subset_25,switch_nodes = c(67),xcross = c(0.999), niter = 1000, cores = 12, model = "poisson_tree", early_growth_model_on = 1.0)
# 
# # Create the plot
# p25 <- ggtree(res_PJBU_25$ultratree, color = "grey") +
#   theme_tree2() +
#   geom_tiplab()+
#   scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 9), breaks = 0:9) +
#   # translucent grey box
#   annotate(
#     "rect",
#     xmin = res_PJBU_25$upper_node_lims$lb95, xmax = res_PJBU_25$upper_node_lims$ub95,
#     ymin = -Inf, ymax = Inf,
#     fill = "grey", alpha = 0.2
#   ) +
#   # dotted vertical lines at x=0 and x=4.1
#   geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
#   geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
#   geom_vline(xintercept = 7.5, color = "red", linetype = "dotted") +
#   geom_vline(xintercept = res_PJBU_25$upper_node_lims$mean, color = "black", linetype = "dotted") +
#   theme(
#     axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
#     axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
#   )
# 
# print(p25)
# 
# 
# tips_to_remove <- c("PJBUGTDABC25","PJBUGTDABC35")
# 
# phylo_tree_subset_25_35 <- drop.tip(phylo_tree, tips_to_remove )
# 
# res_PJBU_25_35=fit_tree(tree=phylo_tree_subset_25_35,switch_nodes = c(66),xcross = c(0.999), niter = 1000, cores = 12, model = "poisson_tree", early_growth_model_on = 1.0)
# 
# # Create the plot
# p25_35<- ggtree(res_PJBU_25_35$ultratree, color = "grey") +
#   theme_tree2() +
#   geom_tiplab()+
#   scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 9), breaks = 0:9) +
#   # translucent grey box
#   annotate(
#     "rect",
#     xmin = res_PJBU_25_35$upper_node_lims$lb95, xmax = res_PJBU_25_35$upper_node_lims$ub95,
#     ymin = -Inf, ymax = Inf,
#     fill = "grey", alpha = 0.2
#   ) +
#   # dotted vertical lines at x=0 and x=4.1
#   geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
#   geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
#   geom_vline(xintercept = 7.5, color = "red", linetype = "dotted") +
#   geom_vline(xintercept = res_PJBU_25_35$upper_node_lims$mean, color = "black", linetype = "dotted") +
#   theme(
#     axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
#     axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
#   )
# 
# print(p25_35)
# 
# # Save as PDF
# ggsave(
#   filename = "/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/Markus_Rtreefit/PJBU_ultrametric_tree_sbs1_blood_node67_0999_01_subset.pdf",
#   plot = p,
#   width = 8,
#   height = 6
# )



# Calculate mutation rates for malignant and normal B-cells

PJBU_malignant_cells <- c( "PJBUGTDABC9","PJBUGTDABC43","PJBUGTDABC54","PJBUGTDABC26",
                           "PJBUGTDABC37","PJBUGTDABC24","PJBUGTDABC11","PJBUGTDABC8","PJBUGTDABC6",
                           "PJBUGTDABC19","PJBUGTDABC42","PJBUGTDABC50","PJBUGTDABC13","PJBUGTDABC29",
                           "PJBUGTDABC30","PJBUGTDABC33","PJBUGTDABC32","PJBUGTDABC51","PJBUGTDABC56",
                           "PJBUGTDABC3","PJBUGTDABC40","PJBUGTDABC27","PJBUGTDABC55","PJBUGTDABC16",
                           "PJBUGTDABC2","PJBUGTDABC31","PJBUGTDABC44","PJBUGTDABC38","PJBUGTDABC39",
                           "PJBUGTDABC21","PJBUGTDABC58","PJBUGTDABC53","PJBUGTDABC49","PJBUGTDABC52",
                           "PJBUGTDABC1","PJBUGTDABC18","PJBUGTDABC45","PJBUGTDABC17","PJBUGTDABC41",
                           "PJBUGTDABC23","PJBUGTDABC35","PJBUGTDABC25") 

PJBU_normal_cells <- c( "PJBUGTDBBC69","PJBUGTDBBC74","PJBUGTDBBC67","PJBUGTDBBC64","PJBUGTDBBC72",
                        "PJBUGTDBBC68","PJBUGTDBBC59","PJBUGTDBBC65","PJBUGTDBBC76",
                        "PJBUGTDBBC71","PJBUGTDBBC81","PJBUGTDBBC70","PJBUGTDBBC63","PJBUGTDBBC80",
                        "PJBUGTDBBC79")

age_PJBU <- 7.5

# Long-form expansion of `samples`, handling "," or "|" and trimming spaces
x_long <- x %>%
  mutate(row_id = row_number()) %>%
  separate_rows(samples, sep = "\\s*[|,]\\s*") %>%
  mutate(samples = trimws(samples)) %>%
  filter(samples != "")

# Keep only the malignant samples of interest
x_long_mal <- x_long %>%
  filter(samples %in% PJBU_malignant_cells)

# Sum corr_branch_lengths across ALL rows where each sample appears
per_sample_totals <- x_long_mal %>%
  group_by(samples) %>%
  summarise(
    total_corr_branch_length = sum(corr_branch_lengths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_corr_branch_length))

# MRCA length

MRCA_length <- x %>%
  filter(branch_id == "m2") %>%
  pull(corr_branch_lengths)

# Get mutations after MRCA

per_sample_totals$MRCA <- MRCA_length

# Remove MRCA mutation load per sample
per_sample_totals <- per_sample_totals %>%
  mutate(total_corr_branch_length_wo_MRCA = total_corr_branch_length - MRCA)

# Compute latency time
latency_time_PJBU <- age_PJBU - res_PJBU_73$upper_node_lims$mean

# Compute mutation rate *per malignant sample*
per_sample_totals <- per_sample_totals %>%
  mutate(mutation_rate_post_expansion = total_corr_branch_length_wo_MRCA / latency_time_PJBU)

# Get mean and SE
mean_mutation_rate_post_expansion_PJBU <- mean(per_sample_totals$mutation_rate_post_expansion, na.rm = TRUE)
se_mutation_rate_post_expansion_PJBU <- sd(per_sample_totals$mutation_rate_post_expansion, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals$mutation_rate_post_expansion)))

mean_mutation_rate_post_expansion_PJBU
se_mutation_rate_post_expansion_PJBU

# Keep only the normal samples 
x_long_norm <- x_long %>%
  filter(samples %in% PJBU_normal_cells)

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
  mutate(mutation_rate = total_corr_branch_length / age_PJBU)

# Get mean and SE
mean_mutation_rate_normal_PJBU <- mean(per_sample_totals_norm$mutation_rate, na.rm = TRUE)
se_mutation_rate_normal_PJBU <- sd(per_sample_totals_norm$mutation_rate, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals_norm$mutation_rate)))

mean_mutation_rate_normal_PJBU
se_mutation_rate_normal_PJBU

# Mutation rate of pre-expansion cell
mutation_rate_pre_expansion_PJBU <- MRCA_length/res_PJBU_73$upper_node_lims$mean
mutation_rate_pre_expansion_PJBU