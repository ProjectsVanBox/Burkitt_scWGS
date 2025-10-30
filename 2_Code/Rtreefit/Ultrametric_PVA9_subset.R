################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to generate ultrametric trees for Donor PVA9
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

tree = readRDS("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PVA9/Filtered_samples_v3/CPW_01/TreeObject0.1.RDS")
#vcf = VariantAnnotation::readVcf("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PVA9/Filtered_samples_v3/PVA9.vep.sub.SNV.autosomal.noBulksNoBlacklist.vcf")

# prepare tree

tree = prepare_tree(tree)

# plot bare tree

plot_gg_tree_base(tree)
plot_gg_tree_base(tree)+ geom_nodelab()
plot_gg_tree(tree, add_branch_length = TRUE, add_bootstrap = F, add_tip_label = F,add_title = "PVA9")

# VCF per branch
PVA9_contribution <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Trees/PVA9_per_branch_sbs_contribution.csv", row.names = 1)
PVA9_contribution$Signature <- rownames(PVA9_contribution)

subset_sum <- PVA9_contribution %>%
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
x[65,] #this row contains NAs (node 15-to-15 connection)
x<- x[-65,] #deleting NA-filled row
phylo_tree$edge.length<-unname(x$branch_length)
phylo_tree$edge
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label,  # Must match tree$tip.label
  age = rep(13.6,126)  # Replace these with your actual ages 6.6 +0.9
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
input_df_PVA9 <- input_df[input_df$Novogene_ID == "PVA9" & input_df$ResolveDNA_version %in% c("v1", "v2.0", "v2"),]
low_callable_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/low_callable_loci.csv')
below_curve_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')
filtered_samples <- unique(c(low_callable_df$Sample_name, below_curve_df$Sample_name, bad_baf_df$Sample_name , fail_vaf_df$samplename))  # samples that didn't pass QC
input_df_PVA9_filtered <- input_df_PVA9 %>% filter(!Sample_name %in% filtered_samples)
input_df_PVA9_filtered$Callable_Loci <- as.numeric(input_df_PVA9_filtered$Callable_Loci)
input_df_PVA9_filtered$sensitivity <- input_df_PVA9_filtered$Callable_Loci / max_callable

# Step 2: Define the product sensitivity calculator
calculate_product_sensitivity <- function(samples_str, sensitivity_df) {
  sample_list <- strsplit(samples_str, "\\|")[[1]]
  sensitivities <- sensitivity_df$sensitivity[sensitivity_df$Sample_name %in% sample_list]
  if (length(sensitivities) == 0) return(1)  # fallback to avoid division by zero
  return(1 - prod(1 - sensitivities, na.rm = TRUE))
}

# Step 3: Apply it to your tibble `x`
#x <- x %>%
#  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_PVA9_filtered),
#         corr_branch_lengths = ctg_branch_lengths / product_sensitivity)

x <- x %>%
  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_PVA9_filtered),
         corr_branch_lengths = sbs1_sbsblood_branch_lengths / product_sensitivity)

x <- x %>%
  mutate(corr_branch_lengths = round(corr_branch_lengths))

phylo_tree$edge.length<-unname(x$corr_branch_lengths)
phylo_tree$edge
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label,  # Must match tree$tip.label
  age = rep(13.6, 126)  # Replace these with your actual ages 12.7 +0.9
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

tips_to_remove <- c("PVA9GTDBBC74") #remove one normal B-cell

phylo_tree_subset_74 <- drop.tip(phylo_tree, tips_to_remove )



# nodes_to_remove <- as.vector(unlist(x[which(x$tip.label %in% tips_to_remove),2]))
# nodes_to_remove <- c(nodes_to_remove,as.vector(unlist(x[which(x$branch_id %in% c("l","o")),"node"])))
# x2 <- x[-which(x$tip.label %in% tips_to_remove),]
# x2 <- x2[-which(x2$branch_id %in% c("l","o")),]
# 
# phylo_tree_subset$edge <- x2[,c(1:2)]
# phylo_tree_subset$edge <- phylo_tree_subset$edge[-which(phylo_tree_subset$edge$parent %in% nodes_to_remove),]

res_PVA9_74=fit_tree(tree=phylo_tree_subset_74,switch_nodes = c(72),xcross = c(0.999), niter = 20000, cores = 12, model = "poisson_tree", early_growth_model_on = 1.0)

# Create the plot
p74 <- ggtree(res_PVA9_74$ultratree, color = "grey") +
  theme_tree2() +
  #geom_tiplab()+
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 14), breaks = 0:14) +
  # translucent grey box
  annotate(
    "rect",
    xmin = res_PVA9_74$upper_node_lims$lb95, xmax = res_PVA9_74$upper_node_lims$ub95,
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.2
  ) +
  # dotted vertical lines at x=0 and x=4.1
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 13.6, color = "red", linetype = "dotted") +
  geom_vline(xintercept = res_PVA9_74$upper_node_lims$mean, color = "black", linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
  )

print(p74)

res_PVA9_74$upper_node_lims$mean
res_PVA9_74$upper_node_lims$ub95
res_PVA9_74$upper_node_lims$lb95

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/PVA9_ultrametric_tree_sbs1_blood_node72_0999_wo_PVA9GTDBBC74.pdf",
  plot = p74,
  width = 8,
  height = 6
)

# Calculate mutation rates for malignant and normal B-cells

PVA9_malignant_cells <- c( "PVA9GTDABC55",
                           "PVA9GTDABC57","PVA9GTDABC54","PVA9GTDABC48","PVA9GTDABC21","PVA9GTDABC15","PVA9GTDABC1",
                           "PVA9GTDABC52","PVA9GTDABC43","PVA9GTDABC2","PVA9GTDABC24","PVA9GTDABC34","PVA9GTDABC61",
                           "PVA9GTDABC60","PVA9GTDABC59","PVA9GTDABC37","PVA9GTDABC56","PVA9GTDABC33","PVA9GTDABC53",
                           "PVA9GTDABC49","PVA9GTDABC32","PVA9GTDABC40","PVA9GTDABC42","PVA9GTDABC44","PVA9GTDABC35",
                           "PVA9GTDABC51","PVA9GTDABC45","PVA9GTDABC39","PVA9GTDABC50","PVA9GTDABC36","PVA9GTDABC62",
                           "PVA9GTDABC58","PVA9GTDABC46","PVA9GTDABC31","PVA9GTDABC47","PVA9GTDABC18","PVA9GTDABC25",
                           "PVA9GTDABC4","PVA9GTDABC16","PVA9GTDABC28","PVA9GTDABC27","PVA9GTDABC11","PVA9GTDABC19",
                           "PVA9GTDABC12","PVA9GTDABC13","PVA9GTDABC23","PVA9GTDABC17","PVA9GTDABC14","PVA9GTDABC10",
                           "PVA9GTDABC6","PVA9GTDABC5","PVA9GTDABC3") 

PVA9_normal_cells <- c("PVA9GTDBBC77","PVA9GTDBBC65","PVA9GTDBBC76","PVA9GTDBBC68","PVA9GTDBBC63",
                       "PVA9GTDBBC67","PVA9GTDBBC70","PVA9GTDBBC69","PVA9GTDBBC73","PVA9GTDBBC66")

age_PVA9 <- 13.6

# Long-form expansion of `samples`, handling "," or "|" and trimming spaces
x_long <- x %>%
  mutate(row_id = row_number()) %>%
  separate_rows(samples, sep = "\\s*[|,]\\s*") %>%
  mutate(samples = trimws(samples)) %>%
  filter(samples != "")

# Keep only the malignant samples of interest
x_long_mal <- x_long %>%
  filter(samples %in% PVA9_malignant_cells)

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
  filter(branch_id == "J3") %>%
  pull(corr_branch_lengths)

# Get mutations after MRCA

per_sample_totals$MRCA <- MRCA_length

# Remove MRCA mutation load per sample
per_sample_totals <- per_sample_totals %>%
  mutate(total_corr_branch_length_wo_MRCA = total_corr_branch_length - MRCA)

# Compute latency time
latency_time_PVA9 <- age_PVA9 - res_PVA9_74$upper_node_lims$mean

# Compute mutation rate *per malignant sample*
per_sample_totals <- per_sample_totals %>%
  mutate(mutation_rate_post_expansion = total_corr_branch_length_wo_MRCA / latency_time_PVA9)

# Get mean and SE
mean_mutation_rate_post_expansion_PVA9 <- mean(per_sample_totals$mutation_rate_post_expansion, na.rm = TRUE)
se_mutation_rate_post_expansion_PVA9 <- sd(per_sample_totals$mutation_rate_post_expansion, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals$mutation_rate_post_expansion)))

mean_mutation_rate_post_expansion_PVA9
se_mutation_rate_post_expansion_PVA9

# Keep only the normal samples 
x_long_norm <- x_long %>%
  filter(samples %in% PVA9_normal_cells)

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
  mutate(mutation_rate = total_corr_branch_length / age_PVA9)

# Get mean and SE
mean_mutation_rate_normal_PVA9 <- mean(per_sample_totals_norm$mutation_rate, na.rm = TRUE)
se_mutation_rate_normal_PVA9 <- sd(per_sample_totals_norm$mutation_rate, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals_norm$mutation_rate)))

mean_mutation_rate_normal_PVA9
se_mutation_rate_normal_PVA9

# Mutation rate of pre-expansion cell
mutation_rate_pre_expansion_PVA9 <- MRCA_length/res_PVA9_74$upper_node_lims$mean
mutation_rate_pre_expansion_PVA9
