################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at mutational signatures of single cell and bulk WGS samples
# Author: Alexander Steemers
# Date: June 2025
################################################################################

# Load libraries

library(MutationalPatterns)
library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome)
library(ChIPpeakAnno)
library(ggplot2)
library(NMF)
library(RColorBrewer)
library(tibble)
library(reshape2)
library(grid)
library(purrr)
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

# Load metadata 

input_df <-  read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
diagnostic_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_overview.csv')
low_callable_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv')
below_curve_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')

# Make blacklist

blacklist <- unique(c(below_curve_df$Sample_name,low_callable_df$Sample_name, bad_baf_df$Sample_name, fail_vaf_df$samplename))

# Load single cell vcf files (make sure to get the right PTATO-filtered file, which depends on which PTA kit version was used)

# Define PTATO directory

ptato_dir <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO"

# Which patients to check

folders_to_check <- c("P3G6", "PRN4", "P856", "PIA9", "PVA9", "PJBU")

# List all snv unfiltered and filtered VCF files 

all_filtered_vcfs <- unlist(lapply(folders_to_check, function(subdir) {
  list.files(
    file.path(ptato_dir, subdir),
    pattern = "snvs.*filtered\\.vcf(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )
}))

# Remove old PTATO vcf file names

all_filtered_vcfs <- all_filtered_vcfs[!grepl("old", all_filtered_vcfs, ignore.case = TRUE)]

# These are samples that have been run through PTAv2 filtering so need to remove the first version of PTA filtering 

samples_to_exclude <- c(
  # Original list
  "PB11197-BLASC-BCELLP1B4",
  "PB11197-BLASC-BCELLP1C4",
  "PB11197-BLASC-BCELLP1I4",
  "PB11197-BLASC-BCELLP1J3",
  "PB11197-BLASC-BCELLP1K4",
  "PB11197-BLASC-BCELLP1L3",
  "PB11197-BLASC-BCELLP1O3",
  "PB11197-BLASC-BCELLP1P3",
  "P3G6GDDABC71",
  "PB14458-BLPL-BCELLP4B3",
  "PB14458-BLPL-BCELLP4B5",
  "PB14458-BLPL-BCELLP4C3",
  "PB14458-BLPL-BCELLP4D3",
  "PB14458-BLPL-BCELLP4D5",
  "PB14458-BLPL-BCELLP4E3",
  "PB14458-BLPL-BCELLP4J3",
  "PB14458-BLPL-BCELLP4K3",
  "PB14458-BLPL-BCELLP4K5",
  "PB14458-BLPL-BCELLP4L3",
  "PB14458-BLPL-BCELLP4L5",
  "PB14458-BLPL-BCELLP4M3",
  "P856GDDUBC32",
  "P856GDDUBC33",
  "P856GDDUBC34",
  "P856GDDUBC40",
  "P856GDDUBC41",
  "P856GDDUBC42",
  "P856GDDUBC43",
  "P856GDDUBC44",
  "P856GDDUBC45",
  "PB14458-BLBM-BCELLP2B3",
  "PB14458-BLBM-BCELLP2B4",
  "PB14458-BLBM-BCELLP2C4",
  "PB14458-BLBM-BCELLP2E4",
  "PB14458-BLBM-BCELLP2F2",
  "PB14458-BLBM-BCELLP2F4",
  "PB14458-BLBM-BCELLP2I2",
  "PB14458-BLBM-BCELLP2L3",
  "PB14458-BLBM-BCELLP2L4",
  "PB14458-BLBM-BCELLP2M4",
  "PB14458-BLBM-BCELLP2N2",
  "PB14458-BLBM-BCELLP2N4",
  "P856GDDBBC46",
  "P856GDDBBC48",
  "P856GDDBBC54",
  "P856GDDBBC57",
  "P856GDDBBC58",
  "P856GDDBBC59",
  "P856GDDBBC60",
  "P856GDDBBC61",
  "P856GDDBBC62",
  "P856GDDBBC63",
  "P856GDDBBC64"
)

# Remove files matching the exclusion rule

pattern_exclude <- paste0(samples_to_exclude, collapse = "|")
pattern_exclude <- paste0("(", pattern_exclude, ").*\\.snvs\\.ptato\\.filtered\\.vcf\\.gz$")

all_filtered_vcfs <- all_filtered_vcfs[!grepl(pattern_exclude, all_filtered_vcfs, ignore.case = TRUE)]

# Filter out blacklist samples

single_cell_sample_names <- sub(".*\\.vep_([^/\\.]+).*", "\\1", all_filtered_vcfs)
scWGS_vcf_files_sub <- all_filtered_vcfs[!single_cell_sample_names %in% blacklist]

# Extract names of single cells

single_cell_sample_names_sub <- single_cell_sample_names[!single_cell_sample_names %in% blacklist]

# Load bulk WGS Burkitt VCF files

bulk_vcf_files <- c(
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID104AAO/SMuRF/PMCID104AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID132AAL/SMuRF/PMCID132AAL.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P3G6/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P3G6_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID821AAL/SMuRF/PMCID821AAL.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID967AAP/SMuRF/PMCID967AAP.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID610AAS/SMuRF/PMCID610AAS.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID867AAT/SMuRF/PMCID867AAT.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID491AAS/SMuRF/PMCID491AAS.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID321AAO/SMuRF/PMCID321AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PVA9/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PVA9_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PRN4/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PRN4_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P856/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P856_D.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P856/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P856_R.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID540AAN/SMuRF/PMCID540AAN.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PJBU/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PJBU_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID509AAT/SMuRF/PMCID509AAT.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID137AAO/SMuRF/PMCID137AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID340AAO/SMuRF/PMCID340AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PIA9/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PIA9_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID690AAT/SMuRF/PMCID690AAT.vep.SMuRF.filtered.sorted.vcf.gz"
  
)

bulk_vcf_files_sub <- bulk_vcf_files[!grepl("690AAT", bulk_vcf_files)] # no myc translocation found in this sample

# Extract bulk names

bulk_sample_names_sub <- c("POFO","POB0","P3G6","PHUI","PUIM","PWSE","PPWW","PJU1","PCUU","PVA9","PRN4","P856_D","P856_R", "PEVR","PJBU","PC1A","P2PS","P2RW","PIA9")

# Convert to novogene names


# Combine all of our own vcf files
vcf_files_pmc <- c(scWGS_vcf_files_sub, bulk_vcf_files_sub)

# Combine names for each vcf file to then convert sample names in my_grl

all_samples <- c(single_cell_sample_names_sub, bulk_sample_names_sub)

# Create mutational matrix and prepare for de novo

my_grl <- read_vcfs_as_granges(vcf_files_pmc, all_samples, ref_genome, type = 'snv')

mut_mat_internal <- mut_matrix(vcf_list = get_mut_type(my_grl, 'snv'), ref_genome = ref_genome)

print(colSums(mut_mat_internal))

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

# combine all datasets and make de novo matrix

mut_mat <- cbind(mut_mat_internal, mutc_mat_external)
print(colSums(mut_mat))

denovo_mat <- mut_mat + 0.0001

# Extract signatures using Rank 3 to 10

# 3 ranks
nmf_res3 <- extract_signatures(denovo_mat, rank = 3, nrun = 100, single_core = TRUE)
colnames(nmf_res3$signatures) <- c("Signature A", "Signature B", "Signature C")
rownames(nmf_res3$contribution) <- c("Signature A", "Signature B", "Signature C")
saveRDS(nmf_res3, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/nmf_res3_", date, ".rds"))

# 4 ranks
nmf_res4 <- extract_signatures(denovo_mat, rank = 4, nrun = 100, single_core = TRUE)
colnames(nmf_res4$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D")
rownames(nmf_res4$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D")
saveRDS(nmf_res4, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/nmf_res4_", date, ".rds"))

# 5 ranks
nmf_res5 <- extract_signatures(denovo_mat, rank = 5, nrun = 100, single_core = TRUE)
colnames(nmf_res5$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
rownames(nmf_res5$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E")
saveRDS(nmf_res5, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/nmf_res5_", date, ".rds"))

# 6 ranks
nmf_res6 <- extract_signatures(denovo_mat, rank = 6, nrun = 100, single_core = TRUE)
colnames(nmf_res6$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
rownames(nmf_res6$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F")
saveRDS(nmf_res6, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/nmf_res6_", date, ".rds"))

# 7 ranks
nmf_res7 <- extract_signatures(denovo_mat, rank = 7, nrun = 100, single_core = TRUE)
colnames(nmf_res7$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
rownames(nmf_res7$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G")
saveRDS(nmf_res7, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/nmf_res7_", date, ".rds"))

# 8 ranks
nmf_res8 <- extract_signatures(denovo_mat, rank = 8, nrun = 100, single_core = TRUE)
colnames(nmf_res8$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                   "Signature H")
rownames(nmf_res8$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                     "Signature H")
saveRDS(nmf_res8, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/nmf_res8_", date, ".rds"))

# 9 ranks
nmf_res9 <- extract_signatures(denovo_mat, rank = 9, nrun = 100, single_core = TRUE)
colnames(nmf_res9$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                   "Signature H", "Signature I")
rownames(nmf_res9$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                     "Signature H", "Signature I")
saveRDS(nmf_res9, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/nmf_res9_", date, ".rds"))

# 10 ranks
nmf_res10 <- extract_signatures(denovo_mat, rank = 10, nrun = 100, single_core = TRUE)
colnames(nmf_res10$signatures) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                    "Signature H", "Signature I", "Signature J")
rownames(nmf_res10$contribution) <- c("Signature A", "Signature B", "Signature C", "Signature D", "Signature E", "Signature F", "Signature G", 
                                      "Signature H", "Signature I", "Signature J")
saveRDS(nmf_res10, file = paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/nmf_res10_", date, ".rds"))

# Read RDS files (need to specify date) to avoid rerunning de novo again

for (i in 3:10) {
  assign(paste0("nmf_res", i),
         readRDS(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/nmf_res", i, "_20250821.rds")))
}


cutoff <- 0.85

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

# 3 ranks
nmf_res3 <- rename_nmf_signatures(nmf_res3, signatures, cutoff = cutoff)
colnames(nmf_res3$signatures)
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank3_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res3$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank3_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank4_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res4$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank4_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank5_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res5$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank5_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank6_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res6$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank6_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank7_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res7$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank7_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank8_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res8$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank8_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank9_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res9$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank9_signature_contribution_", date, ".pdf"))
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
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank10_tri_nuc_profiles_", date, ".pdf"))
plot_96_profile(nmf_res10$signatures, condensed = TRUE)
dev.off()
pdf(paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/rank10_signature_contribution_", date, ".pdf"))
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
  
  out_path <- paste0("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/cosine_heatmap_rank", i, ".pdf")
  
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



plot_1 <- plot_original_vs_reconstructed(
  mut_mat_internal[, 1:269, drop = FALSE],
  nmf_res8$reconstructed[, 1:269, drop = FALSE],
  y_intercept = 0.85,
  ylims = c(0, 1)
) +
  theme(
    axis.text.x = element_text(size = 4)  # smaller numbers → smaller text
  )

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Denovo_inhouse_single_cell_Machado_pcawg/Plots/Cosine_similarity_orig_vs_recon_nmf_rank8.pdf",
  plot = plot_1,
  width = 8,
  height = 6
)

# strict refitting per sample (for internal samples)

required_cols <- c("SBS18", "SBS9", "SBSblood", "SBS17b", "SBS1", "SBS2", "SBS7a")

# From the first signatures object
sigs_from_signatures <- signatures[, required_cols, drop = FALSE]

# From nmf_ref8$signatures
extra_cols <- c("SBSA", "SBSB")
sigs_from_nmf <- nmf_res8$signatures[, extra_cols, drop = FALSE]

# Combine
sigs_to_check <- cbind(sigs_from_signatures, sigs_from_nmf)

colnames(sigs_to_check)

# Bootstrapped refit with selected signatures

all_refits <- list()

for (singlesample in colnames(mut_mat_internal)){
  
  # subset the mut matrix
  mat <- matrix(mut_mat[,singlesample])
  rownames(mat) <- rownames(mut_mat)
  colnames(mat) <- singlesample
  
  # make a dummy column to remove later (needs at least two columns to retain the matrix properties somehow)
  
  # run contri boots
  
  contri_boots <- fit_to_signatures_bootstrapped(mat,
                                                 sigs_to_check,
                                                 n_boots = 100,
                                                 max_delta = 0.002,
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


# Merge bootstrap iterations results

# Loop over each element in the list
for (i in seq_along(all_refits)) {
  current_cols <- colnames(all_refits[[i]])
  missing_cols <- setdiff(colnames(sigs_to_check), current_cols)
  
  # Add each missing column with 0s
  for (col in missing_cols) {
    all_refits[[i]][[col]] <- 0
  }
  
  # Optional: Reorder columns to match the required order (if you want)
  all_refits[[i]] <- all_refits[[i]][, colnames(sigs_to_check)]
}
merged_contriboots <- do.call(rbind, all_refits)
rownames(merged_contriboots) <- sub(".*\\.", "", rownames(merged_contriboots))

contri_tidy <- as.data.frame(merged_contriboots) %>%
  rownames_to_column(var = 'sampleID') %>%
  separate(
    col   = 'sampleID',
    into  = c('sample', 'replicate'),
    sep   = "_(?=[^_]+$)",  # split at the last underscore
    extra = "merge",
    fill  = "right"
  )

contri_tidy2 <- contri_tidy[!colnames(contri_tidy) %in% c("replicate")]

# Summarize and calculate the mean across the required columns
df1 <- contri_tidy2 %>%
  group_by(sample) %>%  # Group by 'sample' column (or any other grouping column)
  summarise(across(colnames(sigs_to_check), mean), .groups = 'drop')

# reorder sample names based on mut_mat_internal

order_samples <- colnames(mut_mat_internal)

df1 <- df1[match(order_samples, df1$sample), ]


df_t1 <- t(df1 %>% column_to_rownames('sample'))
p1_1 <- plot_contribution(df_t1[,1:250],
                        coord_flip = FALSE,
                        mode = "relative"
)
p1_1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p2_1 <- plot_contribution(df_t1[,1:250],
                          coord_flip = FALSE,
                          mode = "absolute"
)
p2_1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p3_bulk <- plot_contribution(df_t1[,251:269],
                          coord_flip = FALSE,
                          mode = "relative"
)
p3_bulk + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Plot results

# plot original versus reconstructed (MutationalPatterns method)
colnames(df_t1) <- colnames(mut_mat_internal)

reconstructed_profiles <- sigs_to_check %*% (df_t1)[,colnames(mut_mat_internal)]

orig_reconstructed_NMF_1 <- plot_original_vs_reconstructed(mut_mat_internal, reconstructed_profiles, 
                                                         y_intercept = 0.85) + 
  theme_CHemALL() + 
  ggTextAxisRotate() +
  geom_bar(stat = 'identity', fill = 'lightblue') +
  geom_hline(yintercept=0.85) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text.x = element_text(size = 4))
print(orig_reconstructed_NMF_1)

# signatures = rows, samples = columns
mtx <- as.matrix(df_t1)

# per-sample totals
col_tot <- colSums(mtx, na.rm = TRUE)

# avoid divide-by-zero
col_tot[col_tot == 0] <- NA

# convert to proportions (0–1) per column
mtx_prop <- sweep(mtx, 2, col_tot, "/")

# set threshold to <5% in all genomes
thr <- 0.05

# row-wise max across samples
row_max <- apply(mtx_prop, 1, function(v) max(v, na.rm = TRUE))

# signatures with <10% in *every* sample
low_everywhere <- row_max < thr & is.finite(row_max)
signatures_low_1 <- rownames(mtx_prop)[low_everywhere]

# results
signatures_low_1

plot_data <- plot_original_vs_reconstructed(
  mut_mat_internal,
  reconstructed_profiles,
  y_intercept = 0.85
)$data

# Adjust the column name if needed — guessing 'cosine_similarity'
cos_vals <- plot_data$cos_sim

# Calculate summary stats
median_cos_1 <- median(cos_vals, na.rm = TRUE)
mean_cos_1 <- mean(cos_vals, na.rm = TRUE)
range_cos_1 <- range(cos_vals, na.rm = TRUE)

mean_cos_1
range_cos_1
median_cos_1

# "SBSA" "SBS2"  "SBSC"  all had <5% contribution per genome so will be removed and strict refittting will be repeated

required_cols <- c("SBS18", "SBS9", "SBSblood", "SBS17b", "SBS1", "SBS7a")

# Define palette
pal <- c("#D2BD96", "#0A9086", "#B3B3B3", "#A62639","#1D3557", "#6C5B78")
names(pal) <- c("SBS1","SBS9","SBS17b","SBS18","SBSblood", "SBS7a")

# From the first signatures object
sigs_from_signatures <- signatures[, required_cols, drop = FALSE]

# Bootstrapped refit with selected signatures

all_refits_2 <- list()

for (singlesample in colnames(mut_mat_internal)){
  
  # subset the mut matrix
  mat <- matrix(mut_mat[,singlesample])
  rownames(mat) <- rownames(mut_mat)
  colnames(mat) <- singlesample
  
  # make a dummy column to remove later (needs at least two columns to retain the matrix properties somehow)
  
  # run contri boots
  
  contri_boots <- fit_to_signatures_bootstrapped(mat,
                                                 sigs_from_signatures,
                                                 n_boots = 100,
                                                 max_delta = 0.002,
                                                 method = "strict"
  )
  
  # set signatures without contribution to 0
  contri_boots <- data.frame(contri_boots)
  for(i in range(1:length(colnames(sigs_from_signatures)))){
    signatu <- colnames(sigs_from_signatures)[i]
    if (!signatu %in% colnames(contri_boots)){
      contri_boots[,signatu] <- 0
    }
  }
  # store result in list
  all_refits_2[[singlesample]] <- contri_boots
  
}

# recommend to save the object here
#saveRDS(all_refits, 'Data/contri_boots_persample_allsamples.RDS')

# Merge bootstrap iterations results

# Loop over each element in the list
for (i in seq_along(all_refits_2)) {
  current_cols <- colnames(all_refits_2[[i]])
  missing_cols <- setdiff(required_cols, current_cols)
  
  # Add each missing column with 0s
  for (col in missing_cols) {
    all_refits_2[[i]][[col]] <- 0
  }
  
  # Optional: Reorder columns to match the required order (if you want)
  all_refits_2[[i]] <- all_refits_2[[i]][, required_cols]
}
merged_contriboots <- do.call(rbind, all_refits_2)
rownames(merged_contriboots) <- sub(".*\\.", "", rownames(merged_contriboots))

contri_tidy <- as.data.frame(merged_contriboots) %>%
  rownames_to_column(var = 'sampleID') %>%
  separate(
    col   = 'sampleID',
    into  = c('sample', 'replicate'),
    sep   = "_(?=[^_]+$)",  # split at the last underscore
    extra = "merge",
    fill  = "right"
  )

contri_tidy2 <- contri_tidy[!colnames(contri_tidy) %in% c("replicate")]

# Summarize and calculate the mean across the required columns
df1 <- contri_tidy2 %>%
  group_by(sample) %>%  # Group by 'sample' column (or any other grouping column)
  summarise(across(required_cols, mean), .groups = 'drop')

df1 <- df1[match(order_samples, df1$sample), ]


df_t1 <- t(df1 %>% column_to_rownames('sample'))
p1_1 <- plot_contribution(df_t1[,1:250],
                          coord_flip = FALSE,
                          mode = "relative",
                          palette = pal
)
p1_1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p2_1 <- plot_contribution(df_t1[,1:250],
                          coord_flip = FALSE,
                          mode = "absolute",
                          palette = pal
)
p2_1 + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p3_bulk <- plot_contribution(df_t1[,251:269],
                             coord_flip = FALSE,
                             mode = "relative",
                             palette = pal
)
p3_bulk + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Get % of mutations per sig

# your contribution matrix: rows = signatures, cols = samples
contrib_mat <- as.matrix(df_t1[, 251:269])

# sanity check: which signatures are present?
rownames(contrib_mat)

# cohort totals per signature (sum across samples)
sig_totals <- rowSums(contrib_mat, na.rm = TRUE)

# grand total across all signatures
all_total <- sum(sig_totals)

# percentages (only if the rows exist)
sbs1_pct     <- if ("SBS1" %in% names(sig_totals))     sig_totals["SBS1"]     / all_total * 100 else NA_real_
sbsblood_pct <- if ("SBSblood" %in% names(sig_totals)) sig_totals["SBSblood"] / all_total * 100 else NA_real_
sbs18_pct <- if ("SBS18" %in% names(sig_totals)) sig_totals["SBS18"] / all_total * 100 else NA_real_
sbs9_pct <- if ("SBS9" %in% names(sig_totals)) sig_totals["SBS9"] / all_total * 100 else NA_real_

# Plot results

# plot original versus reconstructed (MutationalPatterns method)
colnames(df_t1) <- colnames(mut_mat_internal)

reconstructed_profiles_2 <- signatures[,required_cols] %*% (df_t1)[,colnames(mut_mat_internal)]

orig_reconstructed_NMF_2 <- plot_original_vs_reconstructed(mut_mat_internal, reconstructed_profiles_2, 
                                                         y_intercept = 0.8) + 
  theme_CHemALL() + 
  ggTextAxisRotate() +
  geom_bar(stat = 'identity', fill = 'lightblue') +
  geom_hline(yintercept=0.8) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
print(orig_reconstructed_NMF_2)

plot_data <- plot_original_vs_reconstructed(
  mut_mat_internal,
  reconstructed_profiles_2,
  y_intercept = 0.85
)$data

# Adjust the column name if needed — guessing 'cosine_similarity'
cos_vals <- plot_data$cos_sim

# Calculate summary stats
median_cos_2 <- median(cos_vals, na.rm = TRUE)
mean_cos_2 <- mean(cos_vals, na.rm = TRUE)
range_cos_2 <- range(cos_vals, na.rm = TRUE)

mean_cos_2
range_cos_2
median_cos_2

# Now I want to add extra signatures to see if it improves (including PTA artefact signatures)

# ============================================================
# BOOTSTRAPPED "ADD-ONE-SIGNATURE" SWEEP (vs your baseline)
# ============================================================


set.seed(123)  # reproducible bootstraps

# --- Helper: run your exact bootstrap refit and return mean exposures per sample ---
boot_refit_mean_exposures <- function(mut_mat_internal, mut_mat, signatures_subset,
                                      n_boots = 100, max_delta = 0.002, method = "strict",
                                      order_samples = NULL) {
  req_cols <- colnames(signatures_subset)
  all_refits <- vector("list", length = ncol(mut_mat_internal))
  names(all_refits) <- colnames(mut_mat_internal)
  
  for (singlesample in colnames(mut_mat_internal)) {
    # subset single-sample matrix (same as your pipeline)
    mat <- matrix(mut_mat[, singlesample])
    rownames(mat) <- rownames(mut_mat)
    colnames(mat) <- singlesample
    
    contri_boots <- fit_to_signatures_bootstrapped(
      mat, signatures_subset,
      n_boots   = n_boots,
      max_delta = max_delta,
      method    = method
    )
    contri_boots <- data.frame(contri_boots)
    
    # ensure all requested signatures exist in columns (fill missing with 0)
    for (signatu in req_cols) {
      if (!signatu %in% colnames(contri_boots)) {
        contri_boots[[signatu]] <- 0
      }
    }
    all_refits[[singlesample]] <- contri_boots
  }
  
  # harmonize columns & rbind
  for (i in seq_along(all_refits)) {
    missing_cols <- setdiff(req_cols, colnames(all_refits[[i]]))
    for (col in missing_cols) all_refits[[i]][[col]] <- 0
    all_refits[[i]] <- all_refits[[i]][, req_cols, drop = FALSE]
  }
  merged_contriboots <- do.call(rbind, all_refits)
  rownames(merged_contriboots) <- sub(".*\\.", "", rownames(merged_contriboots))
  
  contri_tidy <- as.data.frame(merged_contriboots) %>%
    rownames_to_column(var = "sampleID") %>%
    separate(
      col   = "sampleID",
      into  = c("sample", "replicate"),
      sep   = "_(?=[^_]+$)",
      extra = "merge",
      fill  = "right"
    )
  
  contri_tidy2 <- contri_tidy[!colnames(contri_tidy) %in% c("replicate")]
  
  df_mean <- contri_tidy2 %>%
    group_by(sample) %>%
    summarise(across(all_of(req_cols), mean), .groups = "drop")
  
  if (!is.null(order_samples)) {
    df_mean <- df_mean %>% filter(sample %in% order_samples)
    df_mean <- df_mean[match(order_samples, df_mean$sample), , drop = FALSE]
  }
  
  # return matrix: signatures x samples (like your df_t1)
  expos <- t(df_mean %>% column_to_rownames("sample"))
  return(expos)
}

# --- 1) Baseline exposures (use your df_t1 if already computed) and baseline cosine ---
sigs_base <- signatures[, required_cols, drop = FALSE]

if (exists("df_t1")) {
  df_t1_baseline <- df_t1
} else {
  df_t1_baseline <- boot_refit_mean_exposures(
    mut_mat_internal, mut_mat, sigs_base,
    n_boots = 100, max_delta = 0.002, method = "strict",
    order_samples = if (exists("order_samples")) order_samples else NULL
  )
}

if (exists("plot_data") && "cos_sim" %in% colnames(plot_data)) {
  baseline_cos <- setNames(plot_data$cos_sim, plot_data$sample)
} else {
  reconstructed_baseline <- sigs_base %*% df_t1_baseline
  pd_base <- plot_original_vs_reconstructed(
    mut_mat_internal, reconstructed_baseline, y_intercept = 0.85
  )$data
  baseline_cos <- setNames(pd_base$cos_sim, pd_base$sample)
}

# --- 2) Define extras: all signatures except the baseline ones you used ---
all_sigs <- colnames(signatures)
extras   <- setdiff(all_sigs, required_cols)

# If you want to ignore obvious artifacts etc., you can filter `extras` here, e.g.:
# extras <- setdiff(extras, c("PTA_v1","PTA_v2"))

# --- 3) Sweep: add ONE extra signature at a time (bootstrapped), compute Δcos per sample ---
delta_list <- vector("list", length(extras))
names(delta_list) <- extras

for (sx in extras) {
  sigs_now <- signatures[, c(required_cols, sx), drop = FALSE]
  
  df_t1_now <- boot_refit_mean_exposures(
    mut_mat_internal, mut_mat, sigs_now,
    n_boots = 100, max_delta = 0.002, method = "strict",
    order_samples = if (exists("order_samples")) order_samples else NULL
  )
  
  recon_now <- sigs_now %*% df_t1_now
  pd_now <- plot_original_vs_reconstructed(
    mut_mat_internal, recon_now, y_intercept = 0.85
  )$data
  
  cos_now <- setNames(pd_now$cos_sim, pd_now$sample)
  common_samples <- intersect(names(baseline_cos), names(cos_now))
  
  delta_list[[sx]] <- tibble::tibble(
    sample          = common_samples,
    added_signature = sx,
    cos_sim         = as.numeric(cos_now[common_samples]),
    delta_cos       = as.numeric(cos_now[common_samples] - baseline_cos[common_samples])
  )
}

delta_long <- dplyr::bind_rows(delta_list)

# --- 4) Order x-axis by median Δ for readability & plot ---
sig_order <- delta_long %>%
  group_by(added_signature) %>%
  summarise(median_delta = median(delta_cos, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(median_delta)) %>%
  pull(added_signature)

delta_long <- delta_long %>%
  mutate(added_signature = factor(added_signature, levels = sig_order))

p_delta_box <- ggplot(delta_long, aes(x = added_signature, y = delta_cos)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.alpha = 0.25) +
  labs(
    x = "Extra signature added (bootstrapped; tested one-at-a-time vs baseline)",
    y = "Δ cosine similarity per sample",
    title = "Per-sample Δ cosine similarity by added signature",
    subtitle = "Δ = cos(baseline + 1 extra) − cos(baseline only); exposures are bootstrap means"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p_delta_box)

# --- 5) Optional: summaries you may want to inspect/save ---
sig_summary <- delta_long %>%
  group_by(added_signature) %>%
  summarise(
    mean_delta   = mean(delta_cos, na.rm = TRUE),
    median_delta = median(delta_cos, na.rm = TRUE),
    p95_delta    = quantile(delta_cos, 0.95, na.rm = TRUE),
    n_samples    = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(median_delta))


p_delta_dots <- ggplot(delta_long, aes(x = added_signature, y = delta_cos)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +   # one dot per sample
  labs(
    x = "Extra signature added (bootstrapped, one-at-a-time)",
    y = "Δ cosine similarity per sample",
    title = "Per-sample Δ cosine similarity by extra signature",
    subtitle = "Each dot = one sample (Δ = cos with extra − cos baseline)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p_delta_dots)

p_delta_dots <- ggplot(delta_long, aes(x = added_signature, y = delta_cos)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +   # one dot per sample
  labs(
    x = "Extra signature added (bootstrapped, one-at-a-time)",
    y = "Δ cosine similarity per sample",
    title = "Per-sample Δ cosine similarity by extra signature",
    subtitle = "Each dot = one sample (Δ = cos with extra − cos baseline)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = element_blank()
  )

print(p_delta_dots)

# ============================================================
# BOOTSTRAPPED "REMOVE-ONE-SIGNATURE" SWEEP (vs full baseline)
# ============================================================

# Assumes you still have:
# - required_cols (your baseline)
# - signatures, mut_mat_internal, mut_mat
# - boot_refit_mean_exposures(...) helper from earlier
# - baseline_cos (cosine with full baseline)
# - delta_long from the add-one sweep (columns: sample, added_signature, cos_sim, delta_cos)

# 1) Remove-one sweep
rem_list <- vector("list", length(required_cols))
names(rem_list) <- required_cols

for (sx in required_cols) {
  sigs_now <- signatures[, setdiff(required_cols, sx), drop = FALSE]
  
  df_t1_now <- boot_refit_mean_exposures(
    mut_mat_internal, mut_mat, sigs_now,
    n_boots = 100, max_delta = 0.002, method = "strict",
    order_samples = if (exists("order_samples")) order_samples else NULL
  )
  
  recon_now <- sigs_now %*% df_t1_now
  pd_now <- plot_original_vs_reconstructed(
    mut_mat_internal, recon_now, y_intercept = 0.85
  )$data
  
  cos_now <- setNames(pd_now$cos_sim, pd_now$sample)
  common_samples <- intersect(names(baseline_cos), names(cos_now))
  
  rem_list[[sx]] <- tibble::tibble(
    sample    = common_samples,
    action    = "remove",
    signature = sx,
    cos_sim   = as.numeric(cos_now[common_samples]),
    delta_cos = as.numeric(cos_now[common_samples] - baseline_cos[common_samples])
  )
}

remove_long <- dplyr::bind_rows(rem_list)

# 2) Harmonize the add-one table and combine
adds_long <- delta_long %>%
  dplyr::transmute(
    sample,
    action    = "add",
    signature = as.character(added_signature),
    cos_sim,
    delta_cos
  )

all_delta <- dplyr::bind_rows(adds_long, remove_long)

# 3) (Nice) ordering per panel
sig_order_add <- adds_long %>%
  dplyr::group_by(signature) %>%
  dplyr::summarise(med = median(delta_cos, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(med)) %>%
  dplyr::pull(signature)

sig_order_remove <- remove_long %>%
  dplyr::group_by(signature) %>%
  dplyr::summarise(med = median(delta_cos, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(med)) %>%
  dplyr::pull(signature)

all_delta <- all_delta %>%
  dplyr::mutate(
    signature = dplyr::case_when(
      action == "add"    ~ factor(signature, levels = sig_order_add),
      action == "remove" ~ factor(signature, levels = sig_order_remove)
    )
  )

# 4) One figure: dots for every sample’s Δ, split by action
p_add_remove <- ggplot2::ggplot(all_delta, aes(x = signature, y = delta_cos)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  ggplot2::facet_grid(rows = vars(action), scales = "free_x", space = "free_x") +
  ggplot2::labs(
    x = "Signature",
    y = "Δ cosine similarity per sample",
    title = "Per-sample Δ cosine by signature: add-one vs remove-one",
    subtitle = "Δ = cos(with change) − cos(full baseline); bootstrapped mean exposures"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major.x = ggplot2::element_blank()
  )

print(p_add_remove)

all_delta_one_row <- all_delta %>%
  dplyr::mutate(signature_display = ifelse(action == "add",
                                           paste0("+ ", signature),
                                           paste0("− ", signature)))

ggplot2::ggplot(all_delta_one_row,
                aes(x = signature_display, y = delta_cos, shape = action)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::geom_jitter(width = 0.25, alpha = 0.6, size = 1.6) +
  ggplot2::labs(x = "Change to model (signature)",
                y = "Δ cosine per sample",
                title = "Add-one and Remove-one Δ cosine (bootstrapped)") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))

# Calculate mean ± SE per signature & action
mean_delta <- all_delta_one_row %>%
  group_by(signature_display, action) %>%
  summarise(
    mean_delta = mean(delta_cos, na.rm = TRUE),
    se_delta   = sd(delta_cos, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

pd <- position_dodge(width = 0.5)

ggplot2::ggplot(all_delta_one_row,
                aes(x = signature_display, y = delta_cos, shape = action)) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  #ggplot2::geom_jitter(width = 0.25, alpha = 0.6, size = 1.6) +
  # mean points (no inheritance)
  ggplot2::geom_point(
    data = mean_delta,
    aes(x = signature_display, y = mean_delta, shape = action, color = action, group = action),
    size = 3,
    position = pd,
    inherit.aes = FALSE
  ) +
  # error bars (no inheritance)
  ggplot2::geom_errorbar(
    data = mean_delta,
    aes(x = signature_display,
        ymin = mean_delta - se_delta,
        ymax = mean_delta + se_delta,
        color = action, group = action),
    width = 0.3,
    position = pd,
    inherit.aes = FALSE
  ) +
  ggplot2::labs(x = "Change to model (signature)",
                y = "Δ cosine per sample",
                title = "Add-one and Remove-one Δ cosine (bootstrapped)") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))


output_dir <- "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg"

# Save as RDS (keeps types/attributes exactly)
saveRDS(delta_long,  file.path(output_dir, "add_one_delta_long.rds"))
saveRDS(remove_long, file.path(output_dir, "remove_one_delta_long.rds"))
saveRDS(all_delta,   file.path(output_dir, "add_remove_all_delta.rds"))

