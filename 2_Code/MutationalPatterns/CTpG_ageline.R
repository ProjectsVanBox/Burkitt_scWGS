################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at constant mutational process for timing
# Author: Alexander Steemers
# Date: August 2025
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

setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/")

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

#all_unfiltered_vcfs <- unlist(lapply(folders_to_check, function(subdir) {
#  list.files(
#    file.path(ptato_dir, subdir),
#    pattern = "snvs\\.ptato\\.vcf(\\.gz)?$",
#    recursive = TRUE,
#    full.names = TRUE
#  )
#}))

all_filtered_vcfs <- unlist(lapply(folders_to_check, function(subdir) {
  list.files(
    file.path(ptato_dir, subdir),
    pattern = "snvs.*filtered\\.vcf(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )
}))

# Remove old PTATO vcf file names

#all_unfiltered_vcfs <- all_unfiltered_vcfs[!grepl("old", all_unfiltered_vcfs, ignore.case = TRUE)]
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

# Create mutational matrix and prepare for de novo

my_grl <- read_vcfs_as_granges(single_cell_sample_names_sub, all_samples, ref_genome, type = 'snv')
mut_mat_internal <- mut_matrix(vcf_list = get_mut_type(my_grl, 'snv'), ref_genome = ref_genome)

print(colSums(mut_mat_internal))

normal_b_cells_df <- input_df[input_df$Myc_translocation_IGV == "No" & 
                                input_df$ResolveDNA_version %in% c("v1", "v2", "v2.0"), ]
normal_b_cells_df_sub <- normal_b_cells_df[!normal_b_cells_df$Sample_name %in% blacklist, ]

mut_mat_internal_normal_b <- mut_mat_internal[, colnames(mut_mat_internal) %in% normal_b_cells_df_sub$Sample_name]

mm_ctg <- mut_mat_internal[grep("\\[C>T\\]G", rownames(mut_mat_internal)), ]

mut_mat_internal_normal_b_ctg <- mm_ctg[, colnames(mm_ctg) %in% normal_b_cells_df_sub$Sample_name]

# Sum mutations per sample (column)
mut_load_normal_b <- colSums(mut_mat_internal_normal_b)
mut_load_normal_b_ctg <- colSums(mut_mat_internal_normal_b_ctg)

# Make sure Sample_name is character
normal_b_cells_df_sub$Sample_name <- as.character(normal_b_cells_df_sub$Sample_name)

# Subset age for matching samples
df_normal_b <- data.frame(
  Sample = names(mut_load_normal_b),
  Mutations = mut_load_normal_b
)
df_normal_b$Age <- normal_b_cells_df_sub$Age_at_sampling_Y[match(df_normal_b$Sample, normal_b_cells_df_sub$Sample_name)]

df_normal_b_ctg <- data.frame(
  Sample = names(mut_load_normal_b_ctg),
  Mutations = mut_load_normal_b_ctg
)
df_normal_b_ctg$Age <- normal_b_cells_df_sub$Age_at_sampling_Y[match(df_normal_b_ctg$Sample, normal_b_cells_df_sub$Sample_name)]

ggplot(df_normal_b, aes(x = Age, y = Mutations)) +
  geom_point() +
  geom_jitter(width = 0.4, height = 0.4, size = 2) +
  geom_smooth(method = "lm", color = "darkgreen", se = FALSE) +
  scale_y_continuous(limits = c(0, 2500)) +
  labs(title = "Mutational Load vs Age (Normal B Cells - Total)",
       x = "Age", y = "Mutation Load") +
  theme_CHemALL() +
  theme(
    panel.grid.major.y = element_blank(),  # Remove major horizontal lines
    panel.grid.minor.y = element_blank(),  # Remove minor horizontal lines
    axis.ticks.y = element_blank()         # Remove y-axis ticks
  )

model_total <- lm(Mutations ~ Age, data = df_normal_b)
summary(model_total)

ggplot(df_normal_b_ctg, aes(x = Age, y = Mutations)) +
  geom_point() +
  geom_jitter(width = 0.4, height = 0.4, size = 2) +
  geom_smooth(method = "lm", color = "darkgreen", se = FALSE) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(title = "Mutational Load vs Age (Normal B Cell - C>TpG only)",
       x = "Age", y = "Mutation Load") +
  theme_CHemALL() +
  theme(
    panel.grid.major.y = element_blank(),  # Remove major horizontal lines
    panel.grid.minor.y = element_blank(),  # Remove minor horizontal lines
    axis.ticks.y = element_blank()         # Remove y-axis ticks
  )

model_ctg <- lm(Mutations ~ Age, data = df_normal_b_ctg)
summary(model_ctg)

# Extract the mutation context and type
sbs9_peaks <- c(
  "T[T>G]T",
  "T[T>G]A",
  "T[T>G]G",
  "T[T>C]G"
)
# Subset the matrix
mm_sbs9 <- mut_mat_internal[sbs9_peaks, ]
mut_mat_internal_normal_b_sbs9 <- mm_sbs9[, colnames(mm_sbs9) %in% normal_b_cells_df_sub$Sample_name]

mut_load_normal_b_sbs9 <- colSums(mut_mat_internal_normal_b_sbs9)

df_normal_b_sbs9 <- data.frame(
  Sample = names(mut_load_normal_b_sbs9),
  Mutations = mut_load_normal_b_sbs9
)
df_normal_b_sbs9$Age <- normal_b_cells_df_sub$Age_at_sampling_Y[match(df_normal_b_sbs9$Sample, normal_b_cells_df_sub$Sample_name)]


ggplot(df_normal_b_sbs9, aes(x = Age, y = Mutations)) +
  geom_point() +
  geom_jitter(width = 0.4, height = 0.4, size = 2) +
  geom_smooth(method = "lm", color = "darkgreen", se = FALSE) +
  scale_y_continuous(limits = c(0, 300)) +
  labs(title = "Mutational Load vs Age (Normal B Cell - SBS9 peaks only)",
       x = "Age", y = "Mutation Load") +
  theme_CHemALL() +
  theme(
    panel.grid.major.y = element_blank(),  # Remove major horizontal lines
    panel.grid.minor.y = element_blank(),  # Remove minor horizontal lines
    axis.ticks.y = element_blank()         # Remove y-axis ticks
  )

model_sbs9 <- lm(Mutations ~ Age, data = df_normal_b_sbs9)
summary(model_sbs9)
