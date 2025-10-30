################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at mutational signatures of single cell Machado and bulk WGS samples
# Author: Alexander Steemers
# Date: June 2025
################################################################################

# Load libraries

library(MutationalPatterns)
library(BSgenome)
library(ChIPpeakAnno)
library(ggplot2)
library(NMF)
library(RColorBrewer)
library(tibble)
library(reshape2)
library(grid)
library(readxl)
library(stringr)
library(tidyr)
library(dplyr)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# Load functions and colour palettes

mycols_paired <- brewer.pal(12,"Paired")
mycols_dark2 <- brewer.pal(8, "Dark2")
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Set directory

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/")

# Set date

date <- format(Sys.Date(), "%Y%m%d")

# Import mut matrix from Machado et al. paper (https://www.nature.com/articles/s41586-022-05072-7)

mut_mat_Machado = read.table(file="~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Machado_Nature2020/mutcounts_matrix_AX001_KX001_KX002_KX003_TX001_TX002_CB001 (1).txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
colonyinfo_all = read.table(file="~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Machado_Nature2020/colonyinfo_AX001_KX001_KX002_KX003_TX001_TX002_CB001 (1).txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")
rownames(colonyinfo_all) = colonyinfo_all$colony
colonyinfo_all2 = colonyinfo_all[colnames(mut_mat_Machado),]

# Import mut matrix from PCAWG paper (https://www.nature.com/articles/s41586-020-1969-6#Sec14)

mut_mat_pcawg = read.table("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Machado_Nature2020/mutcounts_matrix_pcawgfocal7_mm_Aug2020 (1).txt", stringsAsFactors = T, header=T)
samppcawg = read.table("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Machado_Nature2020/samplesgroups_mutcounts_matrix_pcawgfocal7_mm_Aug2020.txt", stringsAsFactors = T, header=F)
colnames(samppcawg) = c("samp_names","celltype")
samppcawg$group = samppcawg$celltype
samppcawg$Nmut = apply(mut_mat_pcawg, MARGIN=2, FUN=sum)

# exclude those with more than 60K mutations:

keepsamples = which(samppcawg$Nmut < 60000 & samppcawg$celltype %in% c("Lymph-BNHL","Lymph-CLL","mm","Myeloid-AML") )
mut_mat_pcawg = mut_mat_pcawg[,keepsamples]
samppcawg2 = samppcawg[keepsamples,]

# combine two external datasets

mutc_mat_external = cbind(mut_mat_pcawg, mut_mat_Machado)

denovo_mat <- mutc_mat_external + 0.0001

# Extract signatures using Rank 3 to 10

# 3 ranks
nmf_res3 <- extract_signatures(denovo_mat, rank = 3, nrun = 500, single_core = TRUE)
colnames(nmf_res3$signatures) <- c("Signature A", "Signature B", "Signature C")
rownames(nmf_res3$contribution) <- c("Signature A", "Signature B", "Signature C")
saveRDS(nmf_res3, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_Machado_pcawg/nmf_res3_", date, ".rds"))

# 4 ranks
nmf_res4 <- extract_signatures(denovo_mat, rank = 4, nrun = 500, single_core = TRUE)
colnames(nmf_res4$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
rownames(nmf_res4$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")
saveRDS(nmf_res4, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_Machado_pcawg/nmf_res4_", date, ".rds"))

# 5 ranks
nmf_res5 <- extract_signatures(denovo_mat, rank = 5, nrun = 500, single_core = TRUE)
colnames(nmf_res5$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
rownames(nmf_res5$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
saveRDS(nmf_res5, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_Machado_pcawg/nmf_res5_", date, ".rds"))

# 6 ranks
nmf_res6 <- extract_signatures(denovo_mat, rank = 6, nrun = 500, single_core = TRUE)
colnames(nmf_res6$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
rownames(nmf_res6$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
saveRDS(nmf_res6, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_Machado_pcawg/nmf_res6_", date, ".rds"))

# 7 ranks
nmf_res7 <- extract_signatures(denovo_mat, rank = 7, nrun = 500, single_core = TRUE)
colnames(nmf_res7$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
rownames(nmf_res7$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
saveRDS(nmf_res7, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_Machado_pcawg/nmf_res7_", date, ".rds"))

# 8 ranks
nmf_res8 <- extract_signatures(denovo_mat, rank = 8, nrun = 500, single_core = TRUE)
colnames(nmf_res8$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                   "Signature H")
rownames(nmf_res8$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                     "Signature H")
saveRDS(nmf_res8, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_Machado_pcawg/nmf_res8_", date, ".rds"))

# 9 ranks
nmf_res9 <- extract_signatures(denovo_mat, rank = 9, nrun = 500, single_core = TRUE)
colnames(nmf_res9$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                   "Signature H", "Signature I")
rownames(nmf_res9$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                     "Signature H", "Signature I")
saveRDS(nmf_res9, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_Machado_pcawg/nmf_res9_", date, ".rds"))

# 10 ranks
nmf_res10 <- extract_signatures(denovo_mat, rank = 10, nrun = 500, single_core = TRUE)
colnames(nmf_res10$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                    "Signature H", "Signature I", "Signature J")
rownames(nmf_res10$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                      "Signature H", "Signature I", "Signature J")
saveRDS(nmf_res10, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_Machado_pcawg/nmf_res10_", date, ".rds"))

# Read RDS files (need to specify date) to avoid rerunning de novo again

for (i in 3:10) {
  assign(paste0("nmf_res", i),
         readRDS(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_Machado_pcawg/nmf_res", i, "_20250807.rds")))
}

cutoff <- 0.85

for (i in 3:10) {
  obj_name <- paste0("nmf_res", i)
  nmf_res <- get(obj_name)
  
  rownames(nmf_res$signatures) <- pta_v1_sig[, 1][-length(pta_v1_sig[, 1])]
  
  assign(obj_name, nmf_res)
}

for (i in 3:10) {
  obj_name <- paste0("nmf_res", i)
  nmf_res <- get(obj_name)
  
  new_rownames <- pta_v1_sig[, 1][-length(pta_v1_sig[, 1])]
  
  rownames(nmf_res$signatures) <- new_rownames
  rownames(nmf_res$reconstructed) <- new_rownames
  
  assign(obj_name, nmf_res)
}


# 3 ranks
nmf_res3 <- rename_nmf_signatures(nmf_res3, signatures, cutoff = cutoff)
colnames(nmf_res3$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank3_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res3$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank3_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res3$contribution, nmf_res3$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res3$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 4 ranks
nmf_res4 <- rename_nmf_signatures(nmf_res4, signatures, cutoff = cutoff)
colnames(nmf_res4$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank4_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res4$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank4_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res4$contribution, nmf_res4$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res4$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 5 ranks
nmf_res5 <- rename_nmf_signatures(nmf_res5, signatures, cutoff = cutoff)
colnames(nmf_res5$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank5_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res5$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank5_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res5$contribution, nmf_res5$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res5$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 6 ranks
nmf_res6 <- rename_nmf_signatures(nmf_res6, signatures, cutoff = cutoff)
colnames(nmf_res6$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank6_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res6$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank6_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res6$contribution, nmf_res6$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res6$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 7 ranks
nmf_res7 <- rename_nmf_signatures(nmf_res7, signatures, cutoff = cutoff)
colnames(nmf_res7$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank7_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res7$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank7_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res7$contribution, nmf_res7$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res7$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 8 ranks
nmf_res8 <- rename_nmf_signatures(nmf_res8, signatures, cutoff = cutoff)
colnames(nmf_res8$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank8_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res8$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank8_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res8$contribution, nmf_res8$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res8$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 9 ranks
nmf_res9 <- rename_nmf_signatures(nmf_res9, signatures, cutoff = cutoff)
colnames(nmf_res9$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank9_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res9$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank9_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res9$contribution, nmf_res9$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res9$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 10 ranks
nmf_res10 <- rename_nmf_signatures(nmf_res10, signatures, cutoff = cutoff)
colnames(nmf_res10$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank10_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res10$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/rank10_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res10$contribution, nmf_res10$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res10$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

for (i in 3:10) {
  nmf_res <- get(paste0("nmf_res", i))
  cos_sim <- cos_sim_matrix(nmf_res$signature, signatures)
  
  plot_obj <- plot_cosine_heatmap(
    cos_sim,
    cluster_rows = F,
    cluster_cols = F,
    plot_values = TRUE
  )
  
  out_path <- paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_Machado_pcawg/Plots/cosine_heatmap_rank", i, ".pdf")
  
  ggsave(out_path, plot = plot_obj, width = 12, height = 10)
}

# Loop over values from 3 to 10 to use nmf_res3 to nmf_res10
for (i in 3:10) {
  # Dynamically access nmf_res object based on the loop index
  nmf_res_name <- paste0("nmf_res", i)
  current_nmf_res <- get(nmf_res_name)  # Use get() to access the object dynamically
  
  # Extract signatures for the current NMF result
  current_signatures <- current_nmf_res$signature
  
  # Calculate cosine similarity for the current signatures
  cos_sim_samples_signatures <- cos_sim_matrix(current_signatures, signatures)
  
  # Plot the cosine similarity heatmap for the current iteration
  heatmap_plot <- plot_cosine_heatmap(cos_sim_samples_signatures, 
                                      cluster_rows = TRUE, 
                                      cluster_cols = TRUE, 
                                      plot_values = TRUE)
  print(heatmap_plot)
}
