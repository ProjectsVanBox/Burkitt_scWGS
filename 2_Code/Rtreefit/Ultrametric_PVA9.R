################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to generate ultrametric trees for Donor PVA9
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

tree = readRDS("~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/TreeBuilding_Alex/PVA9/Filtered_samples_v3/CPW_02/TreeObject0.2.RDS")
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
  age = rep(13.6,126)  # Replace these with your actual ages 12.7 +0.9
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
input_df_PVA9 <- input_df[input_df$Novogene_ID == "PVA9" & input_df$ResolveDNA_version %in% c("v1", "v2.0", "v2"),]
low_callable_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv')
below_curve_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')
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
  age = rep(13.6,126)  # Replace these with your actual ages 12.7 +0.9
)
phylo_tree$agedf<-agedf
phylo_tree$edge.length[is.na(phylo_tree$edge.length)] <- 0

p <- ggtree(phylo_tree) +
  geom_text2(aes(label = round(branch.length, 2)), hjust = -0.3) +
  theme_tree2()

print(p)

res_PVA9=fit_tree(tree=phylo_tree,switch_nodes = c(73),xcross = c(0.999), niter = 20000, cores = 4, model = "poisson_tree", early_growth_model_on = 1.0)

# Save ultrametric tree

save(res_PVA9, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/res_PVA9.RData")
load("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/res_PVA9.RData")

# Create the plot
p <- ggtree(res_PVA9$ultratree, color = "grey") +
  theme_tree2() +
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 15), breaks = 0:15) +
  # translucent grey box
  annotate(
    "rect",
    xmin = res_PVA9$upper_node_lims$lb95, xmax = res_PVA9$upper_node_lims$ub95,
    ymin = -Inf, ymax = Inf,
    fill = "grey", alpha = 0.2
  ) +
  # dotted vertical lines at x=0 and x=4.1
  geom_vline(xintercept = 0, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted") +
  geom_vline(xintercept = 13.6, color = "red", linetype = "dotted") +
  geom_vline(xintercept = res_PVA9$upper_node_lims$mean, color = "black", linetype = "dotted") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_text(size = 12, angle = 180, vjust = 0.5)
  )

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/PVA9_ultrametric_tree_sbs1_blood_node73_0999.pdf",
  plot = p,
  width = 8,
  height = 6
)

res_PVA9_74=fit_tree(tree=phylo_tree,switch_nodes = c(74),xcross = c(0.999), niter = 20000, cores = 4, model = "poisson_tree", early_growth_model_on = 1.0)

# Save ultrametric tree

save(res_PVA9_74, file = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/res_PVA9_74.RData")
load("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/res_PVA9_74.RData")

# Create the plot
p <- ggtree(res_PVA9_74$ultratree, color = "grey") +
  theme_tree2() +
  scale_x_continuous(name = "Age (Years since conception)", limits = c(0, 15), breaks = 0:15) +
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

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/PVA9_ultrametric_tree_sbs1_blood_node74_0999.pdf",
  plot = p,
  width = 8,
  height = 6
)