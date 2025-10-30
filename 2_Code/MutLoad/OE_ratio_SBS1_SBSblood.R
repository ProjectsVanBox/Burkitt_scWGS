################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: O/E ratio after subsetting for SBS1 and SBSblood
# Author: Alexander Steemers
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
library(rstatix)
library(tidyr)
library(ggpubr)
library(dplyr)
# Mixed models & helpers
suppressPackageStartupMessages({
  library(lme4)        # lmer
  library(lmerTest)    # p-values for lmer
  library(broom)       # for lm; broom.mixed if available
})

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)

# Load functions and colour palettes
mycols_paired <- brewer.pal(12,"Paired")
mycols_dark2  <- brewer.pal(8, "Dark2")
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Set directory
setwd("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/")

# --------------------------- Load metadata & QC -------------------------------
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

# -------------------- Load single-cell VCF files (PTATO filtered) ------------
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

# ------------------------ Build mutational matrix (autosomes) -----------------
my_grl <- read_vcfs_as_granges(scWGS_vcf_files_sub, single_cell_sample_names_sub, ref_genome, type = 'snv')

autosomes <- paste0("chr", 1:22)
my_grl_autosomes <- endoapply(my_grl, function(x) x[seqnames(x) %in% autosomes])

mut_mat_internal <- mut_matrix(vcf_list = get_mut_type(my_grl_autosomes, 'snv'),
                               ref_genome = ref_genome)

# ---------------------- Keep only non-Myc translocated columns ----------------
input_df_filtered <- input_df[input_df$ResolveDNA_version %in% c("v1", "v2.0", "v2") &
                                !input_df$Sample_name %in% blacklist, , drop = FALSE]
cols_to_keep <- input_df_filtered$Sample_name   # keep all allowed; we’ll label Myc status later
mut_mat_filtered <- mut_mat_internal[, intersect(cols_to_keep, colnames(mut_mat_internal)), drop = FALSE]

# ---------------------- Signature fitting ----------------------
all_signatures <- get_known_signatures()
sbsblood <- read.table("~/Downloads/sigfit_cosmic3_bloodsig_Aug2020.txt", sep = "\t", header = TRUE) |> as.matrix()
SBSblood <- as.numeric(sbsblood[, "Signature.Blood"])
signatures <- cbind(SBSblood = SBSblood, all_signatures)

input_df_filtered <- input_df[input_df$ResolveDNA_version %in% c("v1","v2.0","v2") &
                                !input_df$Sample_name %in% blacklist, , drop = FALSE]
cols_to_keep <- input_df_filtered$Sample_name
mut_mat_filtered <- mut_mat_internal[, intersect(cols_to_keep, colnames(mut_mat_internal)), drop = FALSE]

sub_sig_full <- signatures[, intersect(colnames(signatures),
                                       c("SBS1","SBS7a","SBS9","SBS17b","SBS18","SBSblood")), drop = FALSE]

fit_strict <- fit_to_signatures_strict(mut_mat_filtered, sub_sig_full, max_delta = 0.002)
contrib <- fit_strict$fit_res$contribution

# ---------------------- Normalize and merge ----------------------
callable_fraction <- setNames(input_df$Callable_fraction, input_df$Sample_name)
callable_fraction <- callable_fraction[colnames(contrib)]
callable_fraction[is.na(callable_fraction) | callable_fraction == 0] <- 1
contrib_norm <- sweep(contrib, 2, callable_fraction, "/")

# ---------------------- Add metadata ----------------------
meta <- input_df %>%
  transmute(Sample = Sample_name,
            Age = as.numeric(Age_at_sampling_Y),
            Myc = Myc_translocation_IGV)
df <- as.data.frame(t(contrib_norm)) %>%
  tibble::rownames_to_column("Sample") %>%
  left_join(meta, by = "Sample") %>%
  mutate(MycStatus = ifelse(Myc == "Yes", "Malignant cells", "WT cells"))

# ---------------------- Define WT subtypes by SBS9 ----------------------
# Use SBS9 contribution threshold to define +/−
threshold <- 0.02  # relative contribution (tune as needed)
df <- df %>%
  mutate(
    SBS9_fraction = SBS9 / rowSums(dplyr::select(., starts_with("SBS")), na.rm = TRUE),
    WT_subtype = case_when(
      MycStatus == "WT cells" & SBS9_fraction > threshold ~ "WT_SBS9+",
      MycStatus == "WT cells" & SBS9_fraction <= threshold ~ "WT_SBS9−",
      MycStatus == "Malignant cells" ~ "Malignant cells"
    )
  )

# Merge SBS1 + SBSblood
df <- df %>%
  mutate(SBS1_SBSblood = SBS1 + SBSblood) %>%
  filter(is.finite(SBS1_SBSblood), is.finite(Age))

# ---------------------- Fit model on WT_SBS9− only ----------------------
df_wt_minus <- df %>% filter(WT_subtype == "WT_SBS9−")

fit_lm <- lm(SBS1_SBSblood ~ Age, data = df_wt_minus)
summary(fit_lm)

# Predict expected for all cells
df$Expected <- as.numeric(predict(fit_lm, newdata = df))
df$OE <- df$SBS1_SBSblood / pmax(df$Expected, .Machine$double.eps)

# ---------------------- Scatterplot (WT− vs WT+) ----------------------
p_scatter_wt <- ggplot(df %>% filter(WT_subtype %in% c("WT_SBS9−","WT_SBS9+")),
                       aes(x = Age, y = SBS1_SBSblood, color = WT_subtype, shape = WT_subtype)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(
    data = df %>% filter(WT_subtype == "WT_SBS9−"),
    aes(x = Age, y = SBS1_SBSblood),
    method = "lm", se = FALSE, lwd = 0.8, color = "black"
  ) +
  labs(
    x = "Age (years)",
    y = "SBS1 + SBSblood (per genome)",
    title = "WT subtypes: SBS1+SBSblood vs Age"
  ) +
  scale_color_manual(values = c("WT_SBS9−" = "#e7872b", "WT_SBS9+" = "#4378bd")) +
  theme_CHemALL() +
  theme(legend.title = element_blank())

print(p_scatter_wt)

# ---------------------- OE ratio boxplot ----------------------
df$WT_subtype <- factor(df$WT_subtype, levels = c("WT_SBS9−", "WT_SBS9+", "Malignant cells"))

# ---- Compute Wilcoxon tests ----
pairwise_tests <- pairwise.wilcox.test(df$OE, df$WT_subtype, p.adjust.method = "BH")

cat("\n### Pairwise Wilcoxon Test Results ###\n")
print(pairwise_tests)

# ---- Build the plot ----
p_oe <- ggplot(df, aes(x = WT_subtype, y = OE, fill = WT_subtype)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, color = "black") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", linewidth = 0.7) +
  
  labs(
    x = NULL,
    y = "Observed / Expected (WT_SBS9− trained)",
    title = "OE ratio for merged SBS1+SBSblood"
  ) +
  scale_fill_manual(values = c("WT_SBS9−" = "#e7872b",
                               "WT_SBS9+" = "#cc6d20",
                               "Malignant cells" = "#4378bd")) +
  scale_y_continuous(breaks = 0:6, limits = c(0, 6)) + 
  theme_CHemALL() +
  theme(
    legend.position = "none",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  ggpubr::stat_compare_means(
    comparisons = list(
      c("WT_SBS9−", "WT_SBS9+"),
      c("WT_SBS9−", "Malignant cells"),
      c("WT_SBS9+", "Malignant cells")
    ),
    method = "wilcox.test",
    label = "p.format",
    size = 3.5
  )

print(p_oe)

ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Figures/OE_boxplot_WTsplit_SBS1plusSBSblood.pdf",
       p_oe, width = 5.5, height = 4)


# --- 1) Ensure Age and Biopsy_type are in df ---
age_map <- input_df_filtered %>%
  transmute(Sample = Sample_name,
            Age_cells = suppressWarnings(as.numeric(Age_at_sampling_Y)),
            Biopsy_type = as.character(Biopsy_type))

df <- df %>%
  dplyr::select(-Age) %>%
  left_join(age_map, by = "Sample") %>%
  rename(Age = Age_cells) %>%
  filter(is.finite(Age))

# --- 2) Build x-axis groups: WT pooled, malignant split by Age + Biopsy_type ---
df <- df %>%
  mutate(
    Group_plot = dplyr::case_when(
      WT_subtype %in% c("WT_SBS9−","WT_SBS9+") ~ WT_subtype,
      WT_subtype == "Malignant cells" ~ paste0("Malignant: Age ", Age, " (", Biopsy_type, ")")
    )
  )

# Order: WT first, then malignant sorted by Age ascending, then Biopsy_type
mal_order <- df %>%
  filter(grepl("^Malignant:", Group_plot)) %>%
  distinct(Group_plot, Age, Biopsy_type) %>%
  arrange(Age, Biopsy_type) %>%
  pull(Group_plot)

df$Group_plot <- factor(df$Group_plot,
                        levels = c("WT_SBS9−", "WT_SBS9+", mal_order))

# --- 3) Colors: fixed for WT; for malignant, color by Biopsy_type (consistent across ages) ---
# palette for biopsy types
biopsy_levels <- sort(unique(na.omit(df$Biopsy_type[df$WT_subtype == "Malignant cells"])))
biopsy_pal <- setNames(RColorBrewer::brewer.pal(max(3, length(biopsy_levels)), "Set2")[seq_along(biopsy_levels)],
                       biopsy_levels)

# Set blue for malignant, keep existing WT colours
cols <- c(
  "WT_SBS9−" = "#e7872b",
  "WT_SBS9+" = "#cc6d20",
  # All malignant groups = blue
  setNames(rep("#4378bd", sum(grepl("^Malignant", df$Group_plot))),
           unique(df$Group_plot[grepl("^Malignant", df$Group_plot)]))
)



# 1) Count usable observations per group (non-missing OE)
group_n <- df %>%
  group_by(Group_plot) %>%
  summarise(n = sum(is.finite(OE)), .groups = "drop")

# 2) Keep groups with n >= 2 for testing (do NOT change your plotted df)
valid_groups <- group_n %>% filter(n >= 2) %>% pull(Group_plot)

# 3) Run pairwise Wilcoxon only on valid groups, using WT_SBS9+ as reference
stat_test <- df %>%
  filter(Group_plot %in% valid_groups) %>%
  pairwise_wilcox_test(
    OE ~ Group_plot,
    ref.group = "WT_SBS9−",
    p.adjust.method = "BH",
    detailed = TRUE
  )

# 4) Keep exactly the comparisons you asked for:
#    - WT_SBS9+ vs WT_SBS9−
#    - WT_SBS9+ vs each Malignant:*
stat_test <- stat_test %>%
  filter(
    (group1 == "WT_SBS9−" & group2 == "WT_SBS9+") |
      (group1 == "WT_SBS9−" & grepl("^Malignant:", group2))
  ) %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Group_plot")

# --- 4) Plot: single boxplot, WT pooled, malignant split by Age+Biopsy_type ---
p_oe_by_age_biopsy <- ggplot(df, aes(x = Group_plot, y = OE, fill = Group_plot)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.9, color = "black") +
  geom_hline(yintercept = 1, color = "black", linetype = "dotted", linewidth = 0.6) +
  labs(
    x = NULL,
    y = "Observed / Expected (trained on WT SBS9- cells)"  ) +
  scale_fill_manual(values = cols, guide = "none") +scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 1)) + 
  theme_CHemALL() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 12),
    axis.text.y = element_text(size = 12)
    )

print(p_oe_by_age_biopsy)

# 5) Annotate your existing plot (hide non-sig if you like)
p_oe_by_age_biopsy +
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    hide.ns = FALSE,
    tip.length = 0.01,
    step.increase = 0.06,
    inherit.aes = FALSE   # <- key fix
  )

ggsave(
  "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutLoad/Figures/OE_boxplot_WTpooled_Malignant_AgeByBiopsy.pdf",
  p_oe_by_age_biopsy, width = 5, height = 4
)

# Subset
df_malignant <- df %>%
  dplyr::filter(WT_subtype == "Malignant cells",
                is.finite(OE), is.finite(Age))

# (Optional) fit summary
fit_mal <- lm(OE ~ Age, data = df_malignant)
summary(fit_mal)

# Scatter + linear model
p_mal_scatter <- ggplot(df_malignant, aes(x = Age, y = OE)) +
  geom_point(alpha = 0.85, size = 2, color = "#4378bd") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(
    x = "Age (years)",
    y = "Observed / Expected (SBS1 + SBSblood)",
    title = "Malignant cells: O/E vs Age"
  ) +
  theme_CHemALL()

print(p_mal_scatter)
p_mal_scatter + ggpubr::stat_cor(method = "pearson", label.y.npc = "top", label.x.npc = "left")









# build the list of comparisons: WT_SBS9+ vs every other x group
target <- "WT_SBS9+"
lvls <- levels(df$Group_plot)
comparisons_vs_wt9p <- lapply(setdiff(lvls, target), function(g) c(target, g))

# (optional) also print the exact p-values to console
p_table <- do.call(rbind, lapply(comparisons_vs_wt9p, function(cp) {
  sub <- df %>% dplyr::filter(Group_plot %in% cp)
  p <- tryCatch(wilcox.test(OE ~ Group_plot, data = sub)$p.value, error = function(e) NA_real_)
  data.frame(group1 = cp[1], group2 = cp[2], p_value = p)
}))
cat("\n## Wilcoxon p-values: WT_SBS9+ vs each group ##\n"); print(p_table)

# add the brackets + p-values to your existing plot
ymax <- max(df$OE, na.rm = TRUE)
p_oe_by_age_biopsy <- p_oe_by_age_biopsy +
  ggpubr::stat_compare_means(
    comparisons = comparisons_vs_wt9p,
    method = "wilcox.test",
    label = "p.format",
    size = 3,
    step.increase = 0.08
  ) +
  expand_limits(y = ymax * 1.15)  # give some headroom for brackets

print(p_oe_by_age_biopsy)
