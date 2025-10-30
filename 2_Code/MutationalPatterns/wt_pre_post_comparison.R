################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to look at mutational signatures of non-Myc translocated single cells
#              to separate SBS9+ and SBS9−, compare SNV loads, plot age-lines per signature,
#              compute regression stats, and fit mixed-effects models
#              for BL vs WT (SBS9-high/low) with three planned contrasts.
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
library(ggpubr)
library(dplyr)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# Load functions and colour palettes
mycols_paired <- brewer.pal(12,"Paired")
mycols_dark2  <- brewer.pal(8, "Dark2")
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Set directory
setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/")

# Load metadata 
input_df       <- read_excel('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Sample_overview.xlsx') 
diagnostic_df  <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/Bulk_sample_overview.csv')
low_callable_df<- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/low_callable_loci.csv')
below_curve_df <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/below_curve_samples.csv')
bad_baf_df     <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/bad_baf_samples.csv')
fail_vaf_df    <- read.csv('~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/QC/Data/PTA_samples_failVAFcheck.txt')

# Make blacklist
blacklist <- unique(c(below_curve_df$Sample_name,
                      low_callable_df$Sample_name,
                      bad_baf_df$Sample_name,
                      fail_vaf_df$samplename))

# Load single cell VCF files (PTATO filtered)
ptato_dir <- "~/hpc/pmc_vanboxtel/projects/Burkitt/3_Output/PTATO"
folders_to_check <- c("P3G6", "PRN4", "P856", "PIA9", "PVA9", "PJBU")

all_filtered_vcfs <- unlist(lapply(folders_to_check, function(subdir) {
  list.files(
    file.path(ptato_dir, subdir),
    pattern = "snvs.*filtered\\.vcf(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )
}))

# Remove old PTATO file names
all_filtered_vcfs <- all_filtered_vcfs[!grepl("old", all_filtered_vcfs, ignore.case = TRUE)]

# Samples re-run with PTAv2 filtering → drop earlier filtered files
samples_to_exclude <- c(
  "PB11197-BLASC-BCELLP1B4","PB11197-BLASC-BCELLP1C4","PB11197-BLASC-BCELLP1I4",
  "PB11197-BLASC-BCELLP1J3","PB11197-BLASC-BCELLP1K4","PB11197-BLASC-BCELLP1L3",
  "PB11197-BLASC-BCELLP1O3","PB11197-BLASC-BCELLP1P3","P3G6GDDABC71",
  "PB14458-BLPL-BCELLP4B3","PB14458-BLPL-BCELLP4B5","PB14458-BLPL-BCELLP4C3",
  "PB14458-BLPL-BCELLP4D3","PB14458-BLPL-BCELLP4D5","PB14458-BLPL-BCELLP4E3",
  "PB14458-BLPL-BCELLP4J3","PB14458-BLPL-BCELLP4K3","PB14458-BLPL-BCELLP4K5",
  "PB14458-BLPL-BCELLP4L3","PB14458-BLPL-BCELLP4L5","PB14458-BLPL-BCELLP4M3",
  "P856GDDUBC32","P856GDDUBC33","P856GDDUBC34","P856GDDUBC40","P856GDDUBC41",
  "P856GDDUBC42","P856GDDUBC43","P856GDDUBC44","P856GDDUBC45",
  "PB14458-BLBM-BCELLP2B3","PB14458-BLBM-BCELLP2B4","PB14458-BLBM-BCELLP2C4",
  "PB14458-BLBM-BCELLP2E4","PB14458-BLBM-BCELLP2F2","PB14458-BLBM-BCELLP2F4",
  "PB14458-BLBM-BCELLP2I2","PB14458-BLBM-BCELLP2L3","PB14458-BLBM-BCELLP2L4",
  "PB14458-BLBM-BCELLP2M4","PB14458-BLBM-BCELLP2N2","PB14458-BLBM-BCELLP2N4",
  "P856GDDBBC46","P856GDDBBC48","P856GDDBBC54","P856GDDBBC57","P856GDDBBC58",
  "P856GDDBBC59","P856GDDBBC60","P856GDDBBC61","P856GDDBBC62","P856GDDBBC63","P856GDDBBC64"
)

pattern_exclude <- paste0("(", paste0(samples_to_exclude, collapse = "|"), ").*\\.snvs\\.ptato\\.filtered\\.vcf\\.gz$")
all_filtered_vcfs <- all_filtered_vcfs[!grepl(pattern_exclude, all_filtered_vcfs, ignore.case = TRUE)]

# Filter out blacklist samples
single_cell_sample_names <- sub(".*\\.vep_([^/\\.]+).*", "\\1", all_filtered_vcfs)
scWGS_vcf_files_sub      <- all_filtered_vcfs[!single_cell_sample_names %in% blacklist]
single_cell_sample_names_sub <- single_cell_sample_names[!single_cell_sample_names %in% blacklist]

# Create mutational matrix
my_grl <- read_vcfs_as_granges(scWGS_vcf_files_sub, single_cell_sample_names_sub, ref_genome, type = 'snv')
mut_mat_internal <- mut_matrix(vcf_list = get_mut_type(my_grl, 'snv'), ref_genome = ref_genome)

# Keep only non-Myc translocated
input_df_filtered <- input_df[input_df$ResolveDNA_version %in% c("v1", "v2.0", "v2") &
                                !input_df$Sample_name %in% blacklist, , drop = FALSE]
cols_to_keep <- input_df_filtered$Sample_name[input_df_filtered$Myc_translocation_IGV == "No"]
mut_mat_filtered <- mut_mat_internal[, cols_to_keep, drop = FALSE]

# Get signatures 
all_signatures <- get_known_signatures()
sbsblood <- read.table("~/Downloads/sigfit_cosmic3_bloodsig_Aug2020.txt", sep = "\t", header = TRUE) |> as.matrix()
SBSblood <- as.numeric(sbsblood[, "Signature.Blood"])

signatures <- cbind(SBSblood, all_signatures)

# Refitting part with chosen signatures
required_cols <- c("SBS1","SBS7a", "SBS9", "SBS17b", "SBS18", "SBSblood")
sub_sig <- signatures[, required_cols]

normal_B_refit <- list()

for (singlesample in colnames(mut_mat_filtered)) {
  mat <- as.matrix(mut_mat_filtered[, singlesample, drop = FALSE])
  colnames(mat) <- singlesample
  
  contri_boots <- fit_to_signatures_bootstrapped(
    mat,
    sub_sig,
    n_boots   = 100,
    max_delta = 0.002,
    method    = "strict"
  )
  contri_boots <- data.frame(contri_boots)
  # Ensure all required columns exist (fill missing with 0)
  for (i in seq_along(required_cols)) {
    signatu <- required_cols[i]
    if (!signatu %in% colnames(contri_boots)) {
      contri_boots[, signatu] <- 0
    }
  }
  # Reorder to required_cols
  contri_boots <- contri_boots[, required_cols]
  normal_B_refit[[singlesample]] <- contri_boots
}

names(normal_B_refit) <- NULL
merged_contriboots <- do.call(rbind, normal_B_refit)

contri_tidy <- as.data.frame(merged_contriboots) %>%
  rownames_to_column(var = 'sampleID') %>%
  separate(
    col   = 'sampleID',
    into  = c('sample', 'replicate'),
    sep   = "_(?=[^_]+$)",  # split at the last underscore
    extra = "merge",
    fill  = "right"
  )

# Drop replicate column safely
contri_tidy2 <- contri_tidy[, !(names(contri_tidy) %in% "replicate"), drop = FALSE]

# Summarise: mean across bootstrap replicates per sample
df1 <- contri_tidy2 %>%
  group_by(sample) %>%
  summarise(across(all_of(required_cols), mean), .groups = 'drop')

df_t1 <- t(df1 %>% column_to_rownames('sample'))

# Quick plots
p1 <- plot_contribution(df_t1[, 1:63], coord_flip = FALSE, mode = "relative") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p2 <- plot_contribution(df_t1[, 1:63], coord_flip = FALSE, mode = "absolute") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Group samples by SBS9 proportion (>10% vs ≤10%)
total <- colSums(df_t1)
sbs9_vals    <- df_t1["SBS9", ]
sbs9_prop    <- sbs9_vals / total

samples_sbs9high <- names(sbs9_prop[sbs9_prop > 0.1])
samples_sbs9low  <- names(sbs9_prop[sbs9_prop <= 0.1])

df_high <- data.frame(Sample_name = samples_sbs9high,
                      SBS9_status = "High",
                      stringsAsFactors = FALSE)

df_low  <- data.frame(Sample_name = samples_sbs9low,
                      SBS9_status = "Low",
                      stringsAsFactors = FALSE)

sbs9_status_df <- rbind(df_high, df_low)

# Export to CSV
write.csv(sbs9_status_df, "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data/Denovo_inhouse_single_cell_Machado_pcawg/Normal_SBS9_status_table.csv", row.names = FALSE)

# SNV load per sample and normalization by callable fraction
snv_load <- colSums(mut_mat_filtered)
callable_fraction <- setNames(input_df$Callable_fraction, input_df$Sample_name)
callable_fraction <- callable_fraction[names(snv_load)]
snv_load_normalized <- snv_load / callable_fraction
names(snv_load_normalized) <- names(snv_load)

# Build plot_df with groups
all_samples <- intersect(names(snv_load_normalized), names(sbs9_prop))
plot_df <- data.frame(
  sample        = all_samples,
  snv_load_norm = snv_load_normalized[all_samples],
  group         = ifelse(all_samples %in% samples_sbs9high, "SBS9-high", "SBS9-low"),
  stringsAsFactors = FALSE
)
plot_df$group <- factor(plot_df$group, levels = c("SBS9-high", "SBS9-low"))

# Box plot + t-test p-value 
p1_box <- ggplot(plot_df, aes(x = group, y = snv_load_norm, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  labs(x = "", y = "SNV load (normalized)", title = "SNV load in samples with SBS9high vs SBS9low") +
  scale_y_continuous(limits = c(0, 3000)) +
  scale_fill_manual(values = c("SBS9-high" = "#e7872b", "SBS9-low" = "#a0501a")) +
  theme_minimal(base_size = 14) +
  theme_CHemALL() +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  stat_compare_means(method = "t.test", label = "p.format",
                     comparisons = list(c("SBS9-high", "SBS9-low")))

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Group_comparisons/SBS_high_vs_low_SNV_load_comparison_WT_only_all_samples.pdf",
  plot = p1_box, width = 8, height = 6
)

# Palette for signatures
pal <- c("#D2BD96","#6C5B78", "#0A9086", "#B3B3B3", "#A62639", "#1D3557")
names(pal) <- c("SBS1","SBS7a", "SBS9","SBS17b","SBS18","SBSblood")

# Correct per signature for callable fraction
callable_fraction <- setNames(input_df$Callable_fraction, input_df$Sample_name)
df_corrected <- sweep(df_t1, 2, callable_fraction[colnames(df_t1)], FUN = "/")

# Add metadata
ages <- setNames(input_df$Age_at_sampling_Y, input_df$Sample_name)
df_long <- as.data.frame(t(df_corrected))
df_long$Sample <- rownames(df_long)
df_long$Age <- ages[df_long$Sample]
df_long$CellType <- ifelse(
  input_df$Myc_translocation_IGV[match(df_long$Sample, input_df$Sample_name)] == "No",
  "WT cells", "Malignant cells"
)

# Long format
df_melt <- melt(df_long, id.vars = c("Sample", "Age", "CellType"),
                variable.name = "Signature", value.name = "SNV_per_genome")

df_melt_clean <- df_melt %>%
  na.omit() %>%
  mutate(Age = as.numeric(Age)) %>%
  filter(is.finite(SNV_per_genome), is.finite(Age))

# Plot age vs SNVs per genome (WT only), faceted by signature
p_age <- ggplot(df_melt_clean %>% filter(CellType == "WT cells"),
                aes(x = Age, y = SNV_per_genome, colour = Signature)) +
  geom_jitter(alpha = 0.8, width = 0.5, height = 0.5) +
  scale_colour_manual(values = pal) +
  xlab("Age") + ylab("SNVs / genome") +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.4) +
  facet_wrap(~Signature, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(colour = 'black', size = 10, face = "plain", hjust = 0),
        legend.position = "none",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Total [C>T]G across 4 contexts (WT only), normalized by callable fraction
ctg_motifs <- c("A[C>T]G", "C[C>T]G", "G[C>T]G", "T[C>T]G")
stopifnot(all(ctg_motifs %in% rownames(mut_mat_filtered)))
ctg_counts_total <- colSums(mut_mat_filtered[ctg_motifs, , drop = FALSE])

samples    <- names(ctg_counts_total)
callable_fraction <- setNames(input_df$Callable_fraction,    input_df$Sample_name)
ages              <- setNames(input_df$Age_at_sampling_Y,    input_df$Sample_name)
myc_status        <- setNames(input_df$Myc_translocation_IGV, input_df$Sample_name)

df_ctg <- data.frame(
  Sample         = samples,
  Age            = as.numeric(ages[samples]),
  CellType       = ifelse(myc_status[samples] == "No", "WT cells", "Malignant cells"),
  CTG_per_genome = as.numeric(ctg_counts_total / callable_fraction[samples]),
  stringsAsFactors = FALSE
) %>%
  filter(CellType == "WT cells") %>%
  filter(is.finite(Age), is.finite(CTG_per_genome))

# Plot [C>T]G vs age (WT)
p_ctg <- ggplot(df_ctg, aes(x = Age, y = CTG_per_genome)) +
  geom_jitter(alpha = 0.8, width = 0.4, height = 0.4, color = "#1D3557") +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.6, color = "#1D3557") +
  coord_cartesian(ylim = c(0, 300)) +
  xlab("Age") + ylab("Total [C>T]G per callable genome") +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text  = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Build faceted dataset (WT signatures + [C>T]G_total) 
df_ctg_long <- df_ctg %>%
  transmute(Sample, Age, CellType, Signature = "[C>T]G_total", SNV_per_genome = CTG_per_genome)

df_plot <- bind_rows(
  df_melt_clean %>% filter(CellType == "WT cells"),
  df_ctg_long
)

pal_plus <- c(pal, "[C>T]G_total" = "#444444")

# Faceted plot with free_y, but clamp [C>T]G_total to 0–300
# We add an invisible layer (geom_blank) ONLY for the "[C>T]G_total" facet to
# force its limits to [0, 300], while keeping the others free.
.blank_age <- suppressWarnings(min(df_plot$Age, na.rm = TRUE))
blank_df <- data.frame(
  Age = rep(.blank_age, 2),
  SNV_per_genome = c(0, 300),
  Signature = "[C>T]G_total"
)

plot <- ggplot(df_plot, aes(x = Age, y = SNV_per_genome, colour = Signature)) +
  # Invisible bounds for the [C>T]G_total facet
  geom_blank(data = blank_df, inherit.aes = FALSE,
             aes(x = Age, y = SNV_per_genome)) +
  geom_jitter(alpha = 0.8, width = 0.4, height = 0.4) +
  geom_smooth(method = "lm", se = FALSE, lwd = 0.4) +
  scale_colour_manual(values = pal_plus) +
  facet_wrap(~ Signature, scales = "free_y", nrow = 1) +
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, by = 4)) +
  xlab("Age") + ylab("SNVs / genome") +
  theme_bw() +
  theme(axis.title = element_text(size = 10),
        axis.text  = element_text(size = 6),
        strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(color = "black", size = 6, hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Group_comparisons/WT_ageline_per_sig_and_ctgpg_all_samples.pdf",
  plot = plot, width = 7, height = 2
)

# Regression stats (WT-only data used in df_plot)
my_sigs <- unique(df_plot$Signature)
regression_stats <- lapply(my_sigs, function(sig) {
  fit <- lm(SNV_per_genome ~ Age, data = filter(df_plot, Signature == sig))
  fit_sum <- summary(fit)
  tibble(
    Signature = sig,
    beta      = fit_sum$coefficients[2, 1],
    p_value   = fit_sum$coefficients[2, 4],
    R2        = fit_sum$r.squared,
    RSE       = fit_sum$sigma
  )
}) %>% bind_rows()
print(regression_stats)

# BL vs WT (SBS9 high/low) mixed-effects EMMs 

suppressPackageStartupMessages({
  library(emmeans)
  library(glmmTMB)
})

# Parameters & order (WT-low → WT-high → BL-Trunk → BL-Intermediate → BL-Private)
blwt_sig_keep      <- c("SBS1", "SBS7a", "SBS17b","SBS18","SBS9","SBSblood")
blwt_keep_cols_bl  <- c("BL_Trunk","BL_Intermediate", "BL_Private")
blwt_keep_cols_all <- c("WT_SBS9_low","WT_SBS9_high","BL_Trunk","BL_Intermediate", "BL_Private")
blwt_threshold     <- 0.10  # 10% to match above

# Read BL per-patient matrices
blwt_dir_path <- "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data"
blwt_files <- list.files(blwt_dir_path, pattern = "^sbs_contr_per_group_.*\\.csv$", full.names = TRUE)
stopifnot(length(blwt_files) > 0)
blwt_get_pid <- function(f) sub("^sbs_contr_per_group_(.*)\\.csv$", "\\1", basename(f))
blwt_patient_mats <- setNames(
  lapply(blwt_files, function(f) read.csv(f, header = TRUE, row.names = 1, check.names = FALSE) |> as.matrix()),
  sapply(blwt_files, blwt_get_pid)
)

# Long BL table + relative within Patient × Group
collapse_group <- function(m, prefix) {
  cols <- grep(paste0("^", prefix), colnames(m), value = TRUE)
  if (length(cols) == 0) {
    return(rep(0, nrow(m)))
  } else {
    return(rowSums(m[, cols, drop = FALSE]))
  }
}

blwt_df_long_bl <- dplyr::bind_rows(lapply(names(blwt_patient_mats), function(pid) {
  m <- blwt_patient_mats[[pid]]
  
  data.frame(
    Patient          = pid,
    Signature        = rownames(m),
    BL_Trunk         = collapse_group(m, "BL-Trunk"),
    BL_Intermediate  = collapse_group(m, "BL-Intermediate"),   # includes -PL/-BM/-LN/etc
    BL_Private       = collapse_group(m, "BL-Private"),
    check.names = FALSE
  )
})) |>
  tidyr::pivot_longer(
    cols = c(BL_Trunk, BL_Intermediate, BL_Private),
    names_to = "Group", values_to = "Absolute"
  ) |>
  dplyr::group_by(Patient, Group) |>
  dplyr::mutate(Relative = Absolute / sum(Absolute)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Group = factor(Group, levels = blwt_keep_cols_bl),  # e.g. c("BL_Trunk","BL_Intermediate","BL_Private")
    Signature = factor(Signature, levels = blwt_sig_keep)
  )

# WT groups from df_t1 + input_df
stopifnot(exists("df_t1"), exists("input_df"))
stopifnot("SBS9" %in% rownames(df_t1))

blwt_tot <- colSums(df_t1)
blwt_sbs9_prop  <- df_t1["SBS9", colnames(df_t1)] / blwt_tot

blwt_myc_col <- if ("Myc_translocation_IGV" %in% names(input_df)) "Myc_translocation_IGV" else {
  if ("Myc_translocated_IGV" %in% names(input_df)) "Myc_translocated_IGV" else NA_character_
}
stopifnot(!is.na(blwt_myc_col))

blwt_wt_status  <- setNames(input_df[[blwt_myc_col]], input_df$Sample_name) # "Yes"/"No"
blwt_is_wt      <- blwt_wt_status[colnames(df_t1)] == "No"
blwt_wt_samples <- names(blwt_sbs9_prop)[blwt_is_wt & is.finite(blwt_sbs9_prop)]
stopifnot(length(blwt_wt_samples) > 0)

blwt_wt_group <- ifelse(blwt_sbs9_prop[blwt_wt_samples] > blwt_threshold, "WT_SBS9_high", "WT_SBS9_low")
blwt_wt_abs   <- df_t1[blwt_sig_keep, blwt_wt_samples, drop = FALSE]

blwt_df_long_wt <- lapply(seq_along(blwt_wt_samples), function(i) {
  smp <- blwt_wt_samples[i]
  grp <- blwt_wt_group[i]
  data.frame(
    Patient   = smp,                      # sample as random-effect cluster
    Signature = rownames(blwt_wt_abs),
    Group     = grp,
    Absolute  = blwt_wt_abs[, smp],
    check.names = FALSE
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::group_by(Patient, Group) |>
  dplyr::mutate(Relative = Absolute / sum(Absolute)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Signature = factor(Signature, levels = blwt_sig_keep),
    Group     = factor(Group, levels = c("WT_SBS9_low","WT_SBS9_high"))
  )

# Combine BL + WT and enforce desired x-order
blwt_df_long_all <- dplyr::bind_rows(blwt_df_long_bl, blwt_df_long_wt) |>
  dplyr::mutate(Group = factor(as.character(Group), levels = blwt_keep_cols_all))

# Beta-GLMM per signature, EMMs + three planned contrasts
blwt_eps <- 1e-6
blwt_df_beta <- blwt_df_long_all |>
  dplyr::mutate(Rel_beta = (Relative * (1 - 2*blwt_eps)) + blwt_eps)

blwt_emm_list <- list(); blwt_p_list <- list()

for (sig in levels(blwt_df_beta$Signature)) {
  d <- blwt_df_beta |> dplyr::filter(Signature == sig) |> droplevels()
  if (dplyr::n_distinct(d$Group) < 2) next
  
  fit <- glmmTMB::glmmTMB(Rel_beta ~ Group + (1|Patient),
                          data = d, family = beta_family(link = "logit"))
  
  emms <- emmeans::emmeans(fit, ~ Group)
  link_sum <- as.data.frame(summary(emms, type = "link", infer = c(TRUE, TRUE)))
  
  eta    <- link_sum$emmean
  seEta  <- link_sum$SE
  p_mean <- plogis(eta)
  gprime <- p_mean * (1 - p_mean)
  seP    <- seEta * gprime
  lower_sym <- p_mean - 1.96 * seP
  upper_sym <- p_mean + 1.96 * seP
  
  blwt_emm_df <- data.frame(
    Group     = link_sum$Group,
    emmean    = pmin(pmax(p_mean, 0), 1),
    lower     = pmin(pmax(lower_sym, 0), 1),
    upper     = pmin(pmax(upper_sym, 0), 1),
    Signature = sig
  )
  blwt_emm_list[[sig]] <- blwt_emm_df
  
  # ---- Planned contrasts (fixed) ----
  levs <- levels(emms)$Group
  wvec <- function(named_weights) {
    z <- setNames(rep(0, length(levs)), levs)
    for (nm in names(named_weights)) if (nm %in% levs) z[nm] <- named_weights[[nm]]
    z
  }
  contr_defs <- list()
  if (all(c("BL_Trunk","BL_Intermediate") %in% levs)) {
    contr_defs[["BL Intermediate - Trunk"]] <- wvec(c("BL_Intermediate"= 1, "BL_Trunk"= -1))
  }
  if (all(c("BL_Trunk","BL_Private") %in% levs)) {
    contr_defs[["BL Private - Trunk"]] <- wvec(c("BL_Private"= 1, "BL_Trunk"= -1))
  }
  if (all(c("BL_Intermediate","BL_Private") %in% levs)) {
    contr_defs[["BL Private - BL Intermediate"]] <- wvec(c("BL_Private"= 1, "BL_Intermediate"= -1))  # FIXED
  }
  if (all(c("WT_SBS9_low","WT_SBS9_high") %in% levs)) {
    contr_defs[["WT High - Low"]] <- wvec(c("WT_SBS9_high"= 1, "WT_SBS9_low"= -1))
  }
  if (all(c("WT_SBS9_high","BL_Trunk") %in% levs)) {
    contr_defs[["WT High - BL Trunk"]] <- wvec(c("WT_SBS9_high"= 1, "BL_Trunk"= -1))
  }
  if (length(contr_defs)) {
    diff_df <- as.data.frame(emmeans::contrast(emms, contr_defs, type = "response"))
    blwt_p_list[[sig]] <- diff_df |>
      dplyr::transmute(Signature = sig, contrast, p_value = p.value)
  }
}

blwt_emm_all <- dplyr::bind_rows(blwt_emm_list) |>
  dplyr::mutate(
    Signature = factor(Signature, levels = blwt_sig_keep),
    Group     = factor(as.character(Group), levels = blwt_keep_cols_all)
  )

blwt_pvals <- dplyr::bind_rows(blwt_p_list) |>
  dplyr::mutate(
    contrast = gsub(" / ", " - ", contrast, fixed = TRUE),
    p_label = paste0("p = ", signif(p_value, 3))
  )

# Annotation x-positions computed from actual factor levels 
lvls <- levels(blwt_emm_all$Group)
xpair <- function(a, b) mean(na.omit(match(c(a, b), lvls)))

x_wt_pair        <- xpair("WT_SBS9_low","WT_SBS9_high")
x_cross_pair     <- xpair("WT_SBS9_high","BL_Trunk")
x_bl_trunk_inter <- xpair("BL_Trunk","BL_Intermediate")
x_bl_trunk_priv  <- xpair("BL_Trunk","BL_Private")
x_bl_inter_priv  <- xpair("BL_Intermediate","BL_Private")

# Build annotation data
annot_wt <- blwt_emm_all |>
  dplyr::filter(Group %in% c("WT_SBS9_low","WT_SBS9_high")) |>
  dplyr::group_by(Signature) |>
  dplyr::summarise(y_pos = pmin(1, max(upper, na.rm = TRUE) * 1.08), .groups = "drop") |>
  dplyr::mutate(x_pos = x_wt_pair, contrast = "WT High - Low")

annot_cross <- blwt_emm_all |>
  dplyr::filter(Group %in% c("WT_SBS9_high","BL_Trunk")) |>
  dplyr::group_by(Signature) |>
  dplyr::summarise(y_pos = pmin(1, max(upper, na.rm = TRUE) * 1.08), .groups = "drop") |>
  dplyr::mutate(x_pos = x_cross_pair, contrast = "WT High - BL Trunk")

annot_bl_trunk_inter <- blwt_emm_all |>
  dplyr::filter(Group %in% c("BL_Trunk","BL_Intermediate")) |>
  dplyr::group_by(Signature) |>
  dplyr::summarise(y_pos = pmin(1, max(upper, na.rm = TRUE) * 1.08), .groups = "drop") |>
  dplyr::mutate(x_pos = x_bl_trunk_inter, contrast = "BL Intermediate - Trunk")

annot_bl_trunk_priv <- blwt_emm_all |>
  dplyr::filter(Group %in% c("BL_Trunk","BL_Private")) |>
  dplyr::group_by(Signature) |>
  dplyr::summarise(y_pos = pmin(1, max(upper, na.rm = TRUE) * 1.08), .groups = "drop") |>
  dplyr::mutate(x_pos = x_bl_trunk_priv, contrast = "BL Private - Trunk")  

annot_bl_inter_priv <- blwt_emm_all |>
  dplyr::filter(Group %in% c("BL_Intermediate","BL_Private")) |>
  dplyr::group_by(Signature) |>
  dplyr::summarise(y_pos = pmin(1, max(upper, na.rm = TRUE) * 1.08), .groups = "drop") |>
  dplyr::mutate(x_pos = x_bl_inter_priv, contrast = "BL Private - BL Intermediate")

blwt_annot <- dplyr::bind_rows(
  annot_wt, annot_cross, annot_bl_trunk_inter, annot_bl_trunk_priv, annot_bl_inter_priv
) |>
  dplyr::left_join(blwt_pvals, by = c("Signature","contrast"))

# Plot
p_blwt <- ggplot(blwt_emm_all, aes(x = Group, y = emmean, color = Signature)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.18) +
  facet_wrap(~ Signature, nrow = 1) +
  scale_color_manual(values = pal, drop = FALSE, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    x = NULL, y = "Relative contribution",
    title = "Estimated proportions by group (WT: SBS9-low/high; BL: Trunk/Intermediate/Private)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold")) +
  geom_text(data = blwt_annot, aes(x = x_pos, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 2)

print(p_blwt)

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Group_comparisons/Relative_sign_contri_per_group_all_samples.pdf",
  plot = p_blwt, width = 10, height = 4
)

