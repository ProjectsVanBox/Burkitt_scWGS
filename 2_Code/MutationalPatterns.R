################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at mutational signatures of single cell and bulk WGS samples
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
signatures <- cbind(PTA_v1,PTA_v2, SBSblood, signatures)

# Load metadata 

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
diagnostic_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_overview.csv')
below_curve_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/VAF/Data/PTA_samples_failVAFcheck.txt')

# Make blacklist

blacklist <- unique(c(below_curve_df$Sample_name, bad_baf_df$Sample_name, fail_vaf_df$samplename))

# Load single cell vcf files (make sure to get the right PTATO-filtered file, which depends on which PTA kit version was used)

vcf_files_P3G6_b1 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/P3G6/batch1/snvs/batch1",
                        pattern = "*/*snvs.ptato.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_P3G6_b1 <- vcf_files_P3G6_b1[str_detect(vcf_files_P3G6_b1, "batch1_P3G6GPDABC26|batch1_P3G6GPDABC27")]  # These are the only two in this batch which were not processed with PTAv2
vcf_files_P3G6_b2 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/P3G6/batch2/snvs/batch2",
                                pattern = "*/*snvs.ptato.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_P3G6_b2 <- vcf_files_P3G6_b2[!str_detect(vcf_files_P3G6_b2, "batch2_PB11197-BLASC-BCELLP1O3|batch2_PB11197-BLASC-BCELLP1P3")]  # These are the only two which were processed with PTAv2
vcf_files_P3G6_ptav2 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/P3G6/PTAv2/",
                                pattern = "*/*snvs.ptatoV2.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_PRN4_b1 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PRN4/batch1/snvs/batch1",
                                pattern = "*/*snvs.ptato.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_PRN4_b2 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PRN4/batch2/snvs/batch2",
                                pattern = "*/*snvs.ptato.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_PRN4_b3 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PRN4/batch3/snvs/batch3",
                                pattern = "*/*snvs.ptato.filtered.vcf.gz$", full.names = TRUE, recursive = TRUE )
vcf_files_P856_ptav2 <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/P856/PTAv2",
                                pattern = "*/*snvs.ptatoV2.filtered.vcf$", full.names = TRUE, recursive = TRUE )

scWGS_vcf_files <- c(vcf_files_P3G6_b1, vcf_files_P3G6_b2, vcf_files_P3G6_ptav2, vcf_files_PRN4_b1, vcf_files_PRN4_b2, vcf_files_PRN4_b3, vcf_files_P856_ptav2)

# Filter out blacklist samples and bulk sample

vcf_files_wo_bulk <- scWGS_vcf_files[!str_detect(scWGS_vcf_files, "PRN4GBDLBC72")] # filter bulk sample
scWGS_vcf_files_sub <- vcf_files_wo_bulk[!grepl(paste(blacklist_df$Sample_name, collapse = "|"), vcf_files_wo_bulk)] # filter blacklisted samples

# Load bulk WGS Burkitt VCF files

vcf_file_PRN4_bulk_lymph_node <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PRN4/batch1/snvs/batch1/PB08410_PRN4.vep_PRN4GBDLBC72/PB08410_PRN4.vep_PRN4GBDLBC72.snvs.ptato.filtered.vcf.gz"
vcf_file_P3G6_bulk_ascites <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB11197_new/intermediate/short_variants/somatic_vcfs/PB11197/231018_HMFreg2086_PB11197.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_PB11197-BLASC-BCELLBULK.SMuRF.filtered.vcf.gz"
vcf_file_P856_bulk_pleura <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB14458/intermediate/short_variants/somatic_vcfs/PB14458/231123_HMFreg2101_PB14458.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_PB14458-BLPL-BCELLBULK.SMuRF.filtered.vcf.gz"
vcf_file_P856_bulk_bone_marrow <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO/PB14458/intermediate/short_variants/somatic_vcfs/PB14458/231123_HMFreg2101_PB14458.vcf.filtered_variants_dbnsfp_CosmicCodingMuts_gonl.snps_indels.r5.liftover.hg38.sorted_PB14458-BLBM-BCELLBULK.SMuRF.filtered.vcf.gz"
bulkWGS_vcf_files <- c(vcf_file_PRN4_bulk_lymph_node, vcf_file_P3G6_bulk_ascites, vcf_file_P856_bulk_pleura, vcf_file_P856_bulk_bone_marrow)

# Load bulk WGS Burkitt VCF files from diagnostics (SMuRF filtered as these are PASS + VAF > 0.15)

bulkWGS_diagnostic_vcf_files <- list.files(path = "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples",
                                           pattern = "*/*/*filtered.sorted.vcf.gz$", full.names = TRUE, recursive = TRUE)
bulkWGS_diagnostic_vcf_files_sub <- bulkWGS_diagnostic_vcf_files[!str_detect(bulkWGS_diagnostic_vcf_files, "PMCID211AAO.vep")] # filter out duplicate

# Combine all of our own vcf files
vcf_files_pmc <- c(scWGS_vcf_files_sub, bulkWGS_vcf_files, bulkWGS_diagnostic_vcf_files_sub)

# Extract names for each vcf file to then convert sample names in my_grl

all_samples <- c(input_df$Sample_name, diagnostic_df$PMC_ID)

sample_lookup <- setNames(all_samples, all_samples)

vcf_sample_map <- sapply(vcf_files_pmc, function(path) {
  match <- all_samples[str_detect(path, stringr::fixed(all_samples))]
  if (length(match) == 1) {
    return(match)
  } else if (length(match) > 1) {
    warning(paste("Multiple matches for path:", path, "| Using first match:", match[1]))
    return(match[1])
  } else {
    return(NA)
  }
}, USE.NAMES = TRUE)

matched_vcf_df <- data.frame(
  Sample_or_PMC_ID = vcf_sample_map,
  File_Path = names(vcf_sample_map),
  stringsAsFactors = FALSE
)

# Create mutational matrix and prepare for de novo

my_grl <- read_vcfs_as_granges(vcf_files_pmc, matched_vcf_df$Sample_or_PMC_ID, ref_genome, type = 'snv')
mut_mat_internal <- mut_matrix(vcf_list = get_mut_type(my_grl, 'snv'), ref_genome = ref_genome)

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

mutc_mat_external = cbind(mut_mat_Machado, mut_mat_pcawg)

# combine all datasets and make de novo matrix

mut_mat <- cbind(mut_mat_internal, mutc_mat_external)
denovo_mat <- mut_mat + 0.0001

# Extract signatures using Rank 3 to 10

# 3 ranks
nmf_res3 <- extract_signatures(denovo_mat, rank = 3, nrun = 500, single_core = TRUE)
colnames(nmf_res3$signatures) <- c("Signature A", "Signature B", "Signature C")
rownames(nmf_res3$contribution) <- c("Signature A", "Signature B", "Signature C")
saveRDS(nmf_res3, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/RDS/nmf_res3_", date, ".rds"))

# 4 ranks
nmf_res4 <- extract_signatures(denovo_mat, rank = 4, nrun = 500, single_core = TRUE)
colnames(nmf_res4$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
rownames(nmf_res4$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")
saveRDS(nmf_res4, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/RDS/nmf_res4_", date, ".rds"))

# 5 ranks
nmf_res5 <- extract_signatures(denovo_mat, rank = 5, nrun = 500, single_core = TRUE)
colnames(nmf_res5$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
rownames(nmf_res5$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
saveRDS(nmf_res5, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/RDS/nmf_res5_", date, ".rds"))

# 6 ranks
nmf_res6 <- extract_signatures(denovo_mat, rank = 6, nrun = 500, single_core = TRUE)
colnames(nmf_res6$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
rownames(nmf_res6$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
saveRDS(nmf_res6, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/RDS/nmf_res6_", date, ".rds"))

# 7 ranks
nmf_res7 <- extract_signatures(denovo_mat, rank = 7, nrun = 500, single_core = TRUE)
colnames(nmf_res7$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
rownames(nmf_res7$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
saveRDS(nmf_res7, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/RDS/nmf_res7_", date, ".rds"))

# 8 ranks
nmf_res8 <- extract_signatures(denovo_mat, rank = 8, nrun = 500, single_core = TRUE)
colnames(nmf_res8$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                   "Signature H")
rownames(nmf_res8$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                     "Signature H")
saveRDS(nmf_res8, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/RDS/nmf_res8_", date, ".rds"))

# 9 ranks
nmf_res9 <- extract_signatures(denovo_mat, rank = 9, nrun = 500, single_core = TRUE)
colnames(nmf_res9$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                   "Signature H", "Signature I")
rownames(nmf_res9$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                     "Signature H", "Signature I")
saveRDS(nmf_res9, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/RDS/nmf_res9_", date, ".rds"))

# 10 ranks
nmf_res10 <- extract_signatures(denovo_mat, rank = 10, nrun = 500, single_core = TRUE)
colnames(nmf_res10$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                    "Signature H", "Signature I", "Signature J")
rownames(nmf_res10$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                      "Signature H", "Signature I", "Signature J")
saveRDS(nmf_res10, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/RDS/nmf_res10_", date, ".rds"))

# Read RDS files (need to specify date) to avoid rerunning de novo again

for (i in 3:10) {
  assign(paste0("nmf_res", i),
         readRDS(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/RDS/nmf_res", i, "_20250624.rds")))
}

# Check whether "de novo" signatures are similar to any of the known signatures

# 3 ranks
nmf_res3 <- rename_nmf_signatures(nmf_res3, signatures, cutoff = 0.85)
colnames(nmf_res3$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank3_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res3$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank3_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res3$contribution, nmf_res3$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res3$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 4 ranks
nmf_res4 <- rename_nmf_signatures(nmf_res4, signatures, cutoff = 0.85)
colnames(nmf_res4$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank4_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res4$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank4_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res4$contribution, nmf_res4$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res4$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 5 ranks
nmf_res5 <- rename_nmf_signatures(nmf_res5, signatures, cutoff = 0.85)
colnames(nmf_res5$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank5_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res5$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank5_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res5$contribution, nmf_res5$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res5$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 6 ranks
nmf_res6 <- rename_nmf_signatures(nmf_res6, signatures, cutoff = 0.85)
colnames(nmf_res6$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank6_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res6$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank6_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res6$contribution, nmf_res6$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res6$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 7 ranks
nmf_res7 <- rename_nmf_signatures(nmf_res7, signatures, cutoff = 0.85)
colnames(nmf_res7$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank7_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res7$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank7_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res7$contribution, nmf_res7$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res7$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 8 ranks
nmf_res8 <- rename_nmf_signatures(nmf_res8, signatures, cutoff = 0.85)
colnames(nmf_res8$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank8_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res8$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank8_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res8$contribution, nmf_res8$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res8$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 9 ranks
nmf_res9 <- rename_nmf_signatures(nmf_res9, signatures, cutoff = 0.85)
colnames(nmf_res9$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank9_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res9$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank9_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res9$contribution, nmf_res9$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res9$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

# 10 ranks
nmf_res10 <- rename_nmf_signatures(nmf_res10, signatures, cutoff = 0.85)
colnames(nmf_res10$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank10_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res10$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Denovo/Plots/rank10_signature_contribution_", date, ".pdf"))
plot_contribution(nmf_res10$contribution, nmf_res10$signature, mode = "relative"
) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
plot_compare_profiles(denovo_mat[, 19],
                      nmf_res10$reconstructed[, 19],
                      profile_names = c("Original", "Reconstructed"),
                      condensed = TRUE
)

cos_sim_samples_signatures <- cos_sim_matrix(nmf_res9$signature, signatures)
plot_cosine_heatmap(cos_sim_samples_signatures, 
                    cluster_rows = TRUE, cluster_cols = TRUE, plot_values = T)

#######################################################################################
################### define input variables here
#######################################################################################

required_cols <- c("SBS9", "SBS1", "SBS8", "SBS7a", "SBS18", "SBS17b", "SBS5", "SBS19", "SBS13")

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
                                                 max_delta = 0.001,
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
#saveRDS(all_refits, 'contri_boots_persample_allsamples.RDS')


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

contri_tidy <- as.data.frame(merged_contriboots) %>% rownames_to_column(var = 'sampleID') %>% separate(col = 'sampleID', into = c('sample','replicate'), sep = '\\.')
contri_tidy2 <- contri_tidy[!colnames(contri_tidy) %in% c("replicate")]

## take average of all boots
#if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")
#library(dplyr)

#contri_tidy2_ann <- merge(contri_tidy2, ann_sample_df_ss, by.x = 'sample', by.y = 'SAMPLE')


##### you can group your samples per sample (df_t1) but also per overarching cell type (df_t) if wanted
df1 <- contri_tidy2 %>% group_by(sample) %>% 
  summarise(across(required_cols, mean), .groups = 'drop')
df_t1 <- t(df1 %>% column_to_rownames('sample'))
p1 <- plot_contribution(df_t1,
                       coord_flip = FALSE,
                       mode = "relative"
)

# Adjust x-axis text angle to vertical
p1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p2 <- plot_contribution(df_t1,
                       coord_flip = FALSE,
                       mode = "absolute"
)

# Adjust x-axis text angle to vertical
p2 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
#######################################################################################
################### Plot results
#######################################################################################

# plot original versus reconstructed (MutationalPatterns method)
reconstructed_profiles <- signatures[,required_cols] %*% (df_t1)[,colnames(mut_mat)]

orig_reconstructed_NMF <- plot_original_vs_reconstructed(mut_mat, reconstructed_profiles, 
                                                         y_intercept = 0.85) + 
  theme_CHemALL() + 
  ggTextAxisRotate() +
  geom_bar(stat = 'identity', fill = 'lightblue') +
  geom_hline(yintercept=0.85)



# Transpose to have samples as rows, signatures as columns
df_prop <- prop.table(as.matrix(df_t1), margin = 2)  # divide each column by its sum

# Convert to data.frame for consistency
df_prop <- as.data.frame(df_prop)

# For each signature (row), calculate in how many samples it is <= 10%
low_contrib_frac <- apply(df_prop, 1, function(x) mean(x <= 0.10))

# Find signatures that are <= 10% in at least 95% of samples
low_contrib_signatures <- names(low_contrib_frac[low_contrib_frac >= 0.90])

# Output result
print("Signatures present at <= 10% in at least 95% of samples:")
print(low_contrib_signatures)






#### plot original vs reconstructed cosim and see which samples don't do well
recon_cosims_mat <- cos_sim_matrix(mut_mat, reconstructed_profiles)
recon_cosims <- as.data.frame(diag(recon_cosims_mat)) %>% rownames_to_column('sample')
colnames(recon_cosims)[2] <- 'cosim_recon'

#recon_cosims_ann <- merge(recon_cosims, ann_sample_df_ss, by.x = 'sample', by.y = 'SAMPLE')

ggplot() +
  geom_bar(data = recon_cosims, aes(x = sample, y = cosim_recon),
           stat = "identity", fill = "skyblue3") + theme_CHemALL() +
  #facet_wrap(~ CellType, scales = 'free', ncol = 4) + theme(axis.title.x=element_blank(),
  #                                                          axis.text.x=element_blank(),
  #                                                          axis.ticks.x=element_blank()) +
  ylim(0, 1) + geom_hline(yintercept = 0.9) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


## see if refit of progenitors is worse than HSCs
#recon_cosims_ann %>% group_by(CellType) %>% summarise(meanRecon = mean(cosim_recon))

#### plot contributions per sample
#g <- plot_contribution(df_t1,
#                       coord_flip = FALSE,
#                       mode = "absolute"
#) + theme_CHemALL() + ggTextAxisRotate() + scale_fill_manual(values = c('#0A9086','#b3b3b3', '#e5e5e5','#A62639'))

#g

#g2 <- plot_contribution(df_t1,
#                        coord_flip = FALSE,
#                        mode = "relative"
#) + theme_CHemALL() + ggTextAxisRotate() + scale_fill_manual(values = c('#0A9086','#b3b3b3', '#e5e5e5','#A62639'))

#g2

#################################################
#### plot contributions per group ###############
#################################################
contri_tidy2_ann <- contri_tidy2 %>%
  mutate(CellType = case_when(
    sample %in% c("PMC25993-LNPTA-BCELL1C3", "PDT7GPDLBC01", "PDT7GPDLBC02", "PDT7GPDLBC03", 
                  "PMC25993-LNPTA-HRS1C5", "PDT7GPDLBC04", "PDT7GPDLBC06", "PDT7GPDLBC05", 
                  "PDT7GPDLHR06", "PZ8JGPDLBC03", "PZ8JGPDLBC04", "PZ8JGPDLBC02", 
                  "PZ8JGPDLBC06", "PZ8JGPDLBC01", "PZ8JGPDLBC05") ~ "Bcell",
    
    sample %in% c("PDT7GPDLHR02", "PB31727-HRSLN-HRSCELLSP2C9", "PDT7GPDLHR09", "PDT7GPDLHR07", 
                  "PDT7GPDLHR13", "PDT7GPDLHR03", "PDT7GPDLHR10", "PDT7GPDLHR08", "PDT7GPDLHR04",
                  "PDT7GPDLHR05", "PDT7GPDLHR11", "PMC25993-LNPTA-HRS2E4", "PMC25993-LNPTA-HRS2F4", 
                  "PDT7GPDLHR14", "PB31727-HRSLN-HRSCELLSP2F9", "PDT7GPDLHR12", "PDT7GPDLHR01", 
                  "PB31727-HRSLN-HRSCELLSP2B10", "PB31727-HRSLN-HRSCELLSP2B9", "PZ8JGPDLHR03", 
                  "PZ8JGPDLHR07", "PZ8JGPDLHR08", "PZ8JGPDLHR06", "PZ8JGPDLHR09", "PZ8JGPDLHR10", 
                  "PZ8JGPDLHR12", "PZ8JGPDLHR05", "PZ8JGPDLHR04") ~ "HRS",
    
    TRUE ~ NA_character_  # for samples not in either group
  ))
# average
df2 <- contri_tidy2_ann %>% group_by(CellType) %>% 
  summarise(across(required_cols, mean), .groups = 'drop')

#transpose data frame
df_t <- t(df2 %>% column_to_rownames('CellType'))
# can define plotting order like this
#df_t <- df_t[,c('HSC','MPP','CMP','GMP')]
#df_t <- df_t[c('SBS18','SBS5','SBS1','HSPC'),]

# save output
#saveRDS(df_t, 'ave_bootstr_refit_CHemALLsigs_toCHemALLl_perDonor_perCelltype.RDS')

g <- plot_contribution(df_t,
                       coord_flip = FALSE,
                       mode = "absolute"
) + theme_CHemALL() + ggTextAxisRotate() + scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#000000"))


g

g2 <- plot_contribution(df_t,
                        coord_flip = FALSE,
                        mode = "relative"
) + theme_CHemALL() + ggTextAxisRotate() + scale_fill_manual(values = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"))

g2

