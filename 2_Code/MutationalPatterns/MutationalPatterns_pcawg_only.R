################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at mutational signatures of bulk pcawg WGS samples
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

mutc_mat_external = cbind(mut_mat_pcawg)

# combine all datasets and make de novo matrix

mut_mat <- cbind(mutc_mat_external)
print(colSums(mut_mat))

denovo_mat <- mutc_mat_external + 0.0001

# Extract signatures using Rank 3 to 10

# 3 ranks
nmf_res3 <- extract_signatures(denovo_mat, rank = 3, nrun = 500, single_core = TRUE)
colnames(nmf_res3$signatures) <- c("Signature A", "Signature B", "Signature C")
rownames(nmf_res3$contribution) <- c("Signature A", "Signature B", "Signature C")
saveRDS(nmf_res3, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_pcawg/nmf_res3_", date, ".rds"))

# 4 ranks
nmf_res4 <- extract_signatures(denovo_mat, rank = 4, nrun = 500, single_core = TRUE)
colnames(nmf_res4$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
rownames(nmf_res4$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")
saveRDS(nmf_res4, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_pcawg/nmf_res4_", date, ".rds"))

# 5 ranks
nmf_res5 <- extract_signatures(denovo_mat, rank = 5, nrun = 500, single_core = TRUE)
colnames(nmf_res5$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
rownames(nmf_res5$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
saveRDS(nmf_res5, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_pcawg/nmf_res5_", date, ".rds"))

# 6 ranks
nmf_res6 <- extract_signatures(denovo_mat, rank = 6, nrun = 500, single_core = TRUE)
colnames(nmf_res6$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
rownames(nmf_res6$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
saveRDS(nmf_res6, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_pcawg/nmf_res6_", date, ".rds"))

# 7 ranks
nmf_res7 <- extract_signatures(denovo_mat, rank = 7, nrun = 500, single_core = TRUE)
colnames(nmf_res7$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
rownames(nmf_res7$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
saveRDS(nmf_res7, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_pcawg/nmf_res7_", date, ".rds"))

# 8 ranks
nmf_res8 <- extract_signatures(denovo_mat, rank = 8, nrun = 500, single_core = TRUE)
colnames(nmf_res8$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                   "Signature H")
rownames(nmf_res8$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                     "Signature H")
saveRDS(nmf_res8, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_pcawg/nmf_res8_", date, ".rds"))

# 9 ranks
nmf_res9 <- extract_signatures(denovo_mat, rank = 9, nrun = 500, single_core = TRUE)
colnames(nmf_res9$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                   "Signature H", "Signature I")
rownames(nmf_res9$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                     "Signature H", "Signature I")
saveRDS(nmf_res9, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_pcawg/nmf_res9_", date, ".rds"))

# 10 ranks
nmf_res10 <- extract_signatures(denovo_mat, rank = 10, nrun = 500, single_core = TRUE)
colnames(nmf_res10$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                    "Signature H", "Signature I", "Signature J")
rownames(nmf_res10$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                      "Signature H", "Signature I", "Signature J")
saveRDS(nmf_res10, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_pcawg/nmf_res10_", date, ".rds"))

# Read RDS files (need to specify date) to avoid rerunning de novo again

for (i in 3:10) {
  assign(paste0("nmf_res", i),
         readRDS(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_pcawg/nmf_res", i, "_20250807.rds")))
}

# Check whether "de novo" signatures are similar to any of the known signatures

# Get signatures 

signatures = get_known_signatures()
pta_v1_sig = read.table("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_lymphoma/1_Input/WGS/PTA_Artefact_Signature.txt", sep = "\t", header = T)
pta_v1_sig = as.matrix(pta_v1_sig)
PTA_v1 <- as.numeric(pta_v1_sig[,"PTA"])
PTA_v1 <- PTA_v1[!is.na(PTA_v1)]
pta_v2_sig = read.table("~/hpc/pmc_vanboxtel/resources/signatures/PTAv2_Artefact_Signature.txt", sep = "\t", header = T) 
pta_v2_sig = as.matrix(pta_v2_sig) 
PTA_v2 <- as.numeric(pta_v2_sig[,"PTAv2"])
sbsblood <- read.table("~/surfdrive/Shared/pmc_vanboxtel/personal/asteemers/sigfit_cosmic3_bloodsig_Aug2020.txt", sep = "\t", header = T)
sbsblood = as.matrix(sbsblood)
SBSblood <- as.numeric(sbsblood[,"Signature.Blood"])
signatures <- cbind(PTA_v1, PTA_v2, SBSblood, signatures)
rownames(signatures) <- pta_v1_sig[, 1][-length(pta_v1_sig[, 1])]


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

cutoff <- 0.85

# 3 ranks
nmf_res3 <- rename_nmf_signatures(nmf_res3, signatures, cutoff = cutoff)
colnames(nmf_res3$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank3_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res3$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank3_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank4_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res4$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank4_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank5_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res5$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank5_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank6_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res6$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank6_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank7_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res7$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank7_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank8_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res8$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank8_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank9_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res9$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank9_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank10_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res10$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/rank10_signature_contribution_", date, ".pdf"))
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
  
  out_path <- paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_pcawg/Plots/cosine_heatmap_rank", i, ".pdf")
  
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

#######################################################################################
################### define input variables here
#######################################################################################

required_cols <- c("SBS18", "SBS9", "SBSblood", "SBS17b", "SBS1", "PTA_v1", "SBS7a")

sigs_to_check <- signatures[,required_cols] 

#######################################################################################
################### Bootstrapped refit with selected signatures
#######################################################################################

all_refits <- list()

for (singlesample in colnames(mut_mat)){
  
  # subset the mut matrix
  mat <- matrix(mut_mat[,singlesample])
  rownames(mat) <- rownames(mut_mat)
  colnames(mat) <- singlesample
  
  # make a dummy column to remove later (needs at least two columns to retain the matrix properties somehow)
  
  # run contri boots
  
  contri_boots <- fit_to_signatures_bootstrapped(mat,
                                                 sigs_to_check,
                                                 n_boots = 100,
                                                 #max_delta = 0.001,
                                                 method = "strict"
  )
  
  # set signatures without contribution to 0
  contri_boots <- data.frame(contri_boots)
  for(i in range(1:length(colnames(sigs_to_check)))){
    signatu <- colnames(sigs_to_check)[i]
    if (!signatu %in% colnames(contri_boots)){
      contri_boots[,signatu] <- 0
    }
  }
  # store result in list
  all_refits[[singlesample]] <- contri_boots
  
}

# recommend to save the object here
#saveRDS(all_refits, 'Data/contri_boots_persample_allsamples.RDS')


#######################################################################################
################### Merge bootstrap iterations results
#######################################################################################

# Loop over each element in the list
for (i in seq_along(all_refits)) {
  current_cols <- colnames(all_refits[[i]])
  missing_cols <- setdiff(required_cols, current_cols)
  
  # Add each missing column with 0s
  for (col in missing_cols) {
    all_refits[[i]][[col]] <- 0
  }
  
  # Optional: Reorder columns to match the required order (if you want)
  all_refits[[i]] <- all_refits[[i]][, required_cols]
}
merged_contriboots <- do.call(rbind, all_refits)

p5g <- plot_bootstrapped_contribution(merged_contriboots,mode = "relative",plot_type = "dotplot")
p5h <- plot_bootstrapped_contribution(merged_contriboots)

contri_tidy <- as.data.frame(merged_contriboots) %>%
  rownames_to_column(var = 'sampleID') %>%
  separate(col = 'sampleID', into = c('sample', 'replicate'), sep = '_', extra = "merge", fill = "right")
contri_tidy2 <- contri_tidy[!colnames(contri_tidy) %in% c("replicate")]

##### you can group your samples per sample (df_t1) but also per overarching cell type (df_t) if wanted
# Save the original order of 'sample' column
original_order <- unique(contri_tidy2$sample)

# Summarize and calculate the mean across the required columns
df1 <- contri_tidy2 %>%
  group_by(sample) %>%  # Group by 'sample' column (or any other grouping column)
  summarise(across(required_cols, mean), .groups = 'drop')

# Reorder based on the original order of 'sample' column

df1 <- df1 %>%
  mutate(sample = original_order) %>%
  arrange(match(sample, original_order))
df1 <- df1 %>%
  mutate(sample = sub("\\..*", "", sample))
df1 <- df1[1:268, ]
df1 <- df1 %>%
  mutate(sample = case_when(
    sample %in% folders_to_check ~ paste0(sample, "_bulk"),
    TRUE ~ sample  # Keep the sample name as is for others
  ))
df_t1 <- t(df1 %>% column_to_rownames('sample'))
p1 <- plot_contribution(df_t1[,1:268],
                        coord_flip = FALSE,
                        mode = "relative"
)

# Adjust x-axis text angle to vertical
p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p2 <- plot_contribution(df_t1[,1:268],
                        coord_flip = FALSE,
                        mode = "absolute"
)

# Adjust x-axis text angle to vertical
p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#######################################################################################
################### Plot results
#######################################################################################

# plot original versus reconstructed (MutationalPatterns method)
reconstructed_profiles <- signatures[,required_cols] %*% (df_t1)[,colnames(mut_mat_internal)]

orig_reconstructed_NMF <- plot_original_vs_reconstructed(mut_mat_internal, reconstructed_profiles, 
                                                         y_intercept = 0.85) + 
  theme_CHemALL() + 
  ggTextAxisRotate() +
  geom_bar(stat = 'identity', fill = 'lightblue') +
  geom_hline(yintercept=0.85) +
  coord_cartesian(ylim = c(0, 1))
print(orig_reconstructed_NMF)