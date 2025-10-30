################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to generate ultrametric trees for Donor PIA9
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

tree = readRDS("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PIA9/Filtered_samples_v4/CPW_01/TreeObject0.1.RDS")

# prepare tree

tree = prepare_tree(tree)

# plot bare tree

plot_gg_tree_base(tree)
plot_gg_tree_base(tree)+ geom_nodelab()
plot_gg_tree(tree, add_branch_length = TRUE, add_bootstrap = F, add_tip_label = F,add_title = "PIA9")

# VCF per branch
PIA9_contribution <- read.csv("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Trees/PIA9_per_branch_sbs_contribution.csv", row.names = 1)
PIA9_contribution$Signature <- rownames(PIA9_contribution)

subset_sum <- PIA9_contribution %>%
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
x[61,] #this row contains NAs (node 15-to-15 connection)
x<- x[-61,] #deleting NA-filled row
phylo_tree$edge.length<-unname(x$branch_length)
phylo_tree$edge
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label,  # Must match tree$tip.label
  age = rep(17.9,118)  # Replace these with your actual ages 6.6 +0.9
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
input_df_PIA9 <- input_df[input_df$Novogene_ID == "PIA9" & input_df$ResolveDNA_version %in% c("v1", "v2.0", "v2"),]
low_callable_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/low_callable_loci.csv')
below_curve_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df <- read.csv('~/hpc/pmc_vanboxtel/projects/Burkitt/surfdrive/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')
filtered_samples <- unique(c(low_callable_df$Sample_name, below_curve_df$Sample_name, bad_baf_df$Sample_name , fail_vaf_df$samplename))  # samples that didn't pass QC
input_df_PIA9_filtered <- input_df_PIA9 %>% filter(!Sample_name %in% filtered_samples)
input_df_PIA9_filtered$Callable_Loci <- as.numeric(input_df_PIA9_filtered$Callable_Loci)
input_df_PIA9_filtered$sensitivity <- input_df_PIA9_filtered$Callable_Loci / max_callable

# Step 2: Define the product sensitivity calculator
calculate_product_sensitivity <- function(samples_str, sensitivity_df) {
  sample_list <- strsplit(samples_str, "\\|")[[1]]
  sensitivities <- sensitivity_df$sensitivity[sensitivity_df$Sample_name %in% sample_list]
  if (length(sensitivities) == 0) return(1)  # fallback to avoid division by zero
  return(1 - prod(1 - sensitivities, na.rm = TRUE))
}

# Step 3: Apply it to your tibble `x`
#x <- x %>%
#  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_PIA9_filtered),
#         corr_branch_lengths = ctg_branch_lengths / product_sensitivity)

x <- x %>%
  mutate(product_sensitivity = sapply(samples, calculate_product_sensitivity, sensitivity_df = input_df_PIA9_filtered),
         corr_branch_lengths = sbs1_sbsblood_branch_lengths / product_sensitivity)

x <- x %>%
  mutate(corr_branch_lengths = round(corr_branch_lengths))

phylo_tree$edge.length<-unname(x$corr_branch_lengths)
phylo_tree$edge
print(phylo_tree)
# Create the age dataframe
agedf <- data.frame(
  tip.label = metadata$tip.label,  # Must match tree$tip.label
  age = rep(17.9, 118)  # Replace these with your actual ages 12.7 +0.9
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

tips_to_remove <- c("PIA9GTDBBC56")

phylo_tree_subset_56 <- drop.tip(phylo_tree, tips_to_remove )



# nodes_to_remove <- as.vector(unlist(x[which(x$tip.label %in% tips_to_remove),2]))
# nodes_to_remove <- c(nodes_to_remove,as.vector(unlist(x[which(x$branch_id %in% c("l","o")),"node"])))
# x2 <- x[-which(x$tip.label %in% tips_to_remove),]
# x2 <- x2[-which(x2$branch_id %in% c("l","o")),]
# 
# phylo_tree_subset$edge <- x2[,c(1:2)]
# phylo_tree_subset$edge <- phylo_tree_subset$edge[-which(phylo_tree_subset$edge$parent %in% nodes_to_remove),]

res_PIA9_56=fit_tree(tree=phylo_tree_subset_56,switch_nodes = c(),xcross = c(), niter = 20000, cores = 12, model = "poisson_tree", early_growth_model_on = 1.0)

# No switch node because we found that the O/E was the same between tumour cells and normal B-cells (both naive and memory)
# Create the plot
p56 <- ggtree(res_PIA9_56$ultratree, color = "grey") +
  theme_tree2() +
  #geom_tiplab()+
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 18), breaks = 0:18) +
  # translucent grey box
  #annotate(
  #  "rect",
  #  xmin = res_PIA9_56$upper_node_lims$lb95, xmax = res_PIA9_56$upper_node_lims$ub95,
  #  ymin = -Inf, ymax = Inf,
  #  fill = "grey", alpha = 0.2
  #) +
  # dotted vertical lines at x=0 and x=4.1
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 17.9, color = "red", linetype = "dotted") +
  #geom_vline(xintercept = res_PIA9_56$upper_node_lims$mean, color = "black", linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
  )

print(p56)

# res_PIA9_23$upper_node_lims$mean
# res_PIA9_23$upper_node_lims$ub95
# res_PIA9_23$upper_node_lims$lb95

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/PIA9_ultrametric_tree_sbs1_blood_wo_PIA9GTDBBC56.pdf",
  plot = p56,
  width = 8,
  height = 6
)

res_PIA9_56_switch=fit_tree(tree=phylo_tree_subset_56,switch_nodes = c(64),xcross = c(0.999), niter = 20000, cores = 12, model = "poisson_tree", early_growth_model_on = 1.0)

# No switch node because we found that the O/E was the same between tumour cells and normal B-cells (both naive and memory)
# Create the plot
p56_switch <- ggtree(res_PIA9_56_switch$ultratree, color = "grey") +
  theme_tree2() +
  #geom_tiplab()+
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 18), breaks = 0:18) +
  # translucent grey box
  annotate(
    "rect",
    xmin = res_PIA9_56_switch$upper_node_lims$lb95, xmax = res_PIA9_56_switch$upper_node_lims$ub95,
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.2
  ) +
  # dotted vertical lines at x=0 and x=4.1
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 17.9, color = "red", linetype = "dotted") +
  geom_vline(xintercept = res_PIA9_56_switch$upper_node_lims$mean, color = "black", linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
  )

print(p56_switch)

res_PIA9_56_switch$upper_node_lims$mean
res_PIA9_56_switch$upper_node_lims$ub95
res_PIA9_56_switch$upper_node_lims$lb95

# Calculate mutation rates for malignant and normal B-cells

PIA9_malignant_cells <- c( "PIA9GTDABC37","PIA9GTDABC7","PIA9GTDABC49","PIA9GTDABC21","PIA9GTDABC33",
                           "PIA9GTDABC19","PIA9GTDABC9","PIA9GTDABC26","PIA9GTDABC44","PIA9GTDABC5",
                           "PIA9GTDABC6","PIA9GTDABC25","PIA9GTDABC48","PIA9GTDABC18","PIA9GTDABC40",
                           "PIA9GTDABC34","PIA9GTDABC50","PIA9GTDABC42","PIA9GTDABC46","PIA9GTDABC47",
                           "PIA9GTDABC43","PIA9GTDABC20","PIA9GTDABC36","PIA9GTDABC24","PIA9GTDABC39",
                           "PIA9GTDABC15","PIA9GTDABC51","PIA9GTDABC23","PIA9GTDABC17","PIA9GTDABC45",
                           "PIA9GTDABC41","PIA9GTDABC27","PIA9GTDABC35","PIA9GTDABC22") 

PIA9_normal_cells <- c("PIA9GTDBBC75","PIA9GTDBBC64","PIA9GTDBBC52","PIA9GTDBBC54","PIA9GTDBBC57",
  "PIA9GTDBBC73","PIA9GTDBBC58","PIA9GTDBBC60","PIA9GTDBBC67","PIA9GTDBBC59",
  "PIA9GTDBBC63","PIA9GTDBBC72","PIA9GTDBBC61","PIA9GTDBBC77","PIA9GTDBBC68",
  "PIA9GTDBBC55","PIA9GTDBBC65","PIA9GTDBBC53","PIA9GTDBBC66",
  "PIA9GTDBBC74","PIA9GTDBBC71","PIA9GTDBBC76")

age_PIA9 <- 17.9

# Long-form expansion of `samples`, handling "," or "|" and trimming spaces
x_long <- x %>%
  mutate(row_id = row_number()) %>%
  separate_rows(samples, sep = "\\s*[|,]\\s*") %>%
  mutate(samples = trimws(samples)) %>%
  filter(samples != "")

# Keep only the malignant samples of interest
x_long_mal <- x_long %>%
  filter(samples %in% PIA9_malignant_cells)

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
  filter(branch_id == "a2") %>%
  pull(corr_branch_lengths)

# Get mutations after MRCA

per_sample_totals$MRCA <- MRCA_length

# Remove MRCA mutation load per sample
per_sample_totals <- per_sample_totals %>%
  mutate(total_corr_branch_length_wo_MRCA = total_corr_branch_length - MRCA)

# Compute latency time
latency_time_PIA9 <- age_PIA9 - 9.797081334

# Compute mutation rate *per malignant sample*
per_sample_totals <- per_sample_totals %>%
  mutate(mutation_rate_post_expansion = total_corr_branch_length_wo_MRCA / latency_time_PIA9)

# Get mean and SE
mean_mutation_rate_post_expansion_PIA9 <- mean(per_sample_totals$mutation_rate_post_expansion, na.rm = TRUE)
se_mutation_rate_post_expansion_PIA9 <- sd(per_sample_totals$mutation_rate_post_expansion, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals$mutation_rate_post_expansion)))

mean_mutation_rate_post_expansion_PIA9
se_mutation_rate_post_expansion_PIA9

# Keep only the normal samples 
x_long_norm <- x_long %>%
  filter(samples %in% PIA9_normal_cells)

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
  mutate(mutation_rate = total_corr_branch_length / age_PIA9)

# Get mean and SE
mean_mutation_rate_normal_PIA9 <- mean(per_sample_totals_norm$mutation_rate, na.rm = TRUE)
se_mutation_rate_normal_PIA9 <- sd(per_sample_totals_norm$mutation_rate, na.rm = TRUE) /
  sqrt(sum(!is.na(per_sample_totals_norm$mutation_rate)))

mean_mutation_rate_normal_PIA9
se_mutation_rate_normal_PIA9

# Mutation rate of pre-expansion cell
mutation_rate_pre_expansion_PIA9 <- MRCA_length/9.797081334
mutation_rate_pre_expansion_PIA9
