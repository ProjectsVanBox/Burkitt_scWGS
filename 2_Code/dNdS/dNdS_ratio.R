################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: Script to compare dN/dS ratios between normal and malignant single-cells
# Author: Alexander Steemers
# Date: August 2025
################################################################################

# Load libraries

library(dplyr)
library(VariantAnnotation)
library(GenomeInfoDb)
library(BiocGenerics)
library(ggplot2)
library(dndscv)
library(purrr)
library(tidyr)
library(tibble)
library(readr)
library(ggrepel)

# Function to generate mutations from vcf lists

vcf_list_to_muts <- function(branch_vcf, patient_id) {
  stopifnot(is.list(branch_vcf), length(branch_vcf) > 0)
  dplyr::bind_rows(lapply(names(branch_vcf), function(bid) {
    v <- branch_vcf[[bid]]
    if (is.null(v) || length(v) == 0) return(NULL)
    rr   <- rowRanges(v)
    alts <- unlist(VariantAnnotation::alt(v))
    data.frame(
      patient  = patient_id,
      branchID = bid,
      chr      = as.character(GenomeInfoDb::seqnames(rr)),
      pos      = BiocGenerics::start(rr),
      ref      = as.character(VariantAnnotation::ref(v)),
      mut      = as.character(alts),
      stringsAsFactors = FALSE
    )
  })) %>% 
    dplyr::filter(chr %in% as.character(1:22)) # autosomes only
}

# Function to apply a patient-specific branch->group mapping

apply_group_map <- function(muts_branch, branch_to_group) {
  muts_branch %>%
    filter(branchID %in% names(branch_to_group)) %>%
    mutate(sampleID = unname(branch_to_group[branchID])) %>%
    dplyr::select(sampleID, chr, pos, ref, mut, patient, branchID)
}

# Extract overall dN/dS ("wall") from a dndscv fit, robust to format

extract_overall <- function(res, g){
  if (is.null(res)) return(NULL)
  tab <- res$globaldnds
  if (is.data.frame(tab) && all(c("name","mle","cilow","cihigh") %in% names(tab))) {
    row <- tab[tab$name == "wall", c("mle","cilow","cihigh")]
    data.frame(Group=g, omega_all=row$mle, lo=row$cilow, hi=row$cihigh)
  } else {
    rn <- rownames(tab); stopifnot("wall" %in% rn)
    data.frame(Group=g, omega_all=tab["wall","mle"], lo=tab["wall","cilow"], hi=tab["wall","cihigh"])
  }
}

# Load RDS lists

branch_vcf_P3G6 <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_vcf_P3G6.rds")
branch_vcf_PVA9 <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_vcf_PVA9.rds")
branch_vcf_PIA9 <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_vcf_PIA9.rds")
branch_vcf_PJBU <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_vcf_PJBU.rds")
branch_vcf_P856 <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_vcf_P856.rds")
branch_vcf_PRN4 <- readRDS("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Data/branch_vcf_PRN4.rds")


# Define branches per sample

# P3G6
group1_P3G6 <- c("Z")
group2_P3G6 <- c("T","R","Y","X","P")
group3_P3G6 <- c("J","I","H","N","M","Q","S","W","V","U")
group4_P3G6 <- c("B","A","D","F","G","a")
branch_to_group_P3G6 <- c(
  setNames(rep("BL-Trunk",        length(group1_P3G6)), group1_P3G6),
  setNames(rep("BL-Intermediate", length(group2_P3G6)), group2_P3G6),
  setNames(rep("BL-Private",      length(group3_P3G6)), group3_P3G6),
  setNames(rep("WT",              length(group4_P3G6)), group4_P3G6)
)

# PVA9
group1_PVA9 <- c("J3")
group2_PVA9 <- c("H3","G3","q2","N","v2","E3")
group3_PVA9 <- c("O2","N2","Q2","S2","H2","G2","J2","L2","V2","F2","C2","B2","A2","a2","Z2","z","e2","u","t","w","r","q","o","n","c","b","e","h","g","k","a","Z","l2","k2","n2","Q","P","O","T","W","V","A3","z2","x2","w2","D3","s2","r2","u2","M","L","I3")
group4_PVA9 <- c("L3","K3","N3","K","G","F","I","T3","S3","B","A")
branch_to_group_PVA9 <- c(
  setNames(rep("BL-Trunk",        length(group1_PVA9)), group1_PVA9),
  setNames(rep("BL-Intermediate", length(group2_PVA9)), group2_PVA9),
  setNames(rep("BL-Private",      length(group3_PVA9)), group3_PVA9),
  setNames(rep("WT",              length(group4_PVA9)), group4_PVA9)
)

# PIA9
group1_PIA9 <- c("a2")
group2_PIA9 <- c("d2","E3","C3","t2","s2","e2","G3","I3","J3")
group3_PIA9 <- c("L","K","N","P","I","H","G","T","W","V","Y","B","A","M3","E","D","b","H3","I2","k2","n2","j2","q2","h2","g2","f2","z2","y2","x2","v2","u2","D3","c2","b2")
group4_PIA9 <- c("s","r","u","u","w","p","o","z","n","m","E2","D2","I","I2","K2","M2","O2","k","R2","j","i","V2","X2","h","c","d")
branch_to_group_PIA9 <- c(
  setNames(rep("BL-Trunk",        length(group1_PIA9)), group1_PIA9),
  setNames(rep("BL-Intermediate", length(group2_PIA9)), group2_PIA9),
  setNames(rep("BL-Private",      length(group3_PIA9)), group3_PIA9),
  setNames(rep("WT",              length(group4_PIA9)), group4_PIA9)
)

# PJBU
group1_PJBU <- c("m2")
group2_PJBU <- c("T","l2","j2","T2","S","Q","L")
group3_PJBU <- c("k","j","m","i","p","h","s","f","e","d","w","c","Z","Y","W","V","G2","F2","E2","D2","K2","M2","O2","C2","B2","A2","U2","W2","b2","a2","Z2","f2","e2","Y2","k2","N","M","P","R","K","J","I")
group4_PJBU <- c("p2","o2","r2","t2","v2","n2","A3","z2","C3","G3","F3","I3","G","F","D","A")
branch_to_group_PJBU <- c(
  setNames(rep("BL-Trunk",        length(group1_PJBU)), group1_PJBU),
  setNames(rep("BL-Intermediate", length(group2_PJBU)), group2_PJBU),
  setNames(rep("BL-Private",      length(group3_PJBU)), group3_PJBU),
  setNames(rep("WT",              length(group4_PJBU)), group4_PJBU)
)

# P856
group1_P856 <- c("o") # BL trunk
group2_P856 <- c("j","p", "r") # BL intermediate 
group3_P856 <- c("c", "b", "e", "h", "g", "J", "I", "L", "H", "O", "G", "R", "U", "T", "W", "Y", "C", "B", "A", "C2", "A2", "F", "u", "t", "w", "q") # BL private 
group4_P856 <- c("l", "k") # WT
branch_to_group_P856 <- c(
  setNames(rep("BL-Trunk",        length(group1_P856)), group1_P856),
  setNames(rep("BL-Intermediate", length(group2_P856)), group2_P856),
  setNames(rep("BL-Private",      length(group3_P856)), group3_P856),
  setNames(rep("WT",              length(group4_P856)), group4_P856)
)

# PRN4
group1_PRN4 <- c("h") # BL trunk 
group2_PRN4 <- c("p", "q", "a", "o", "m") # BL intermediate 
group3_PRN4 <- c("G", "F", "E", "C", "B", "A", "v", "u", "x", "t", "M", "L", "O", "Q", "S", "V", "U", "Z", "Y", "j", "i", "l", "n") # BL private 
group4_PRN4 <- c("g", "b", "c") # WT
branch_to_group_PRN4 <- c(
  setNames(rep("BL-Trunk",        length(group1_PRN4)), group1_PRN4),
  setNames(rep("BL-Intermediate", length(group2_PRN4)), group2_PRN4),
  setNames(rep("BL-Private",      length(group3_PRN4)), group3_PRN4),
  setNames(rep("WT",              length(group4_PRN4)), group4_PRN4)
)

# Build one mutations table across patients

muts_group_all <- bind_rows(
  apply_group_map(vcf_list_to_muts(branch_vcf_P3G6, "P3G6"), branch_to_group_P3G6),
  apply_group_map(vcf_list_to_muts(branch_vcf_PVA9, "PVA9"), branch_to_group_PVA9),
  apply_group_map(vcf_list_to_muts(branch_vcf_PIA9, "PIA9"), branch_to_group_PIA9),
  apply_group_map(vcf_list_to_muts(branch_vcf_PJBU, "PJBU"), branch_to_group_PJBU),
  apply_group_map(vcf_list_to_muts(branch_vcf_P856, "P856"), branch_to_group_P856),
  apply_group_map(vcf_list_to_muts(branch_vcf_PRN4, "PRN4"), branch_to_group_PRN4)
)

# Relabel groups (keep 4 groups)

muts_phase <- muts_group_all %>%
  mutate(sampleID = case_when(
    sampleID == "BL-Trunk"        ~ "BL Pre-expansion",
    sampleID == "BL-Intermediate" ~ "BL-Intermediate",
    sampleID == "BL-Private"      ~ "BL-Private",
    sampleID == "WT"              ~ "WT",
    TRUE ~ sampleID
  ))

# Exclude immunoglobulin loci (hg38): https://pubmed.ncbi.nlm.nih.gov/18432650/#:~:text=The%20human%20immunoglobulins%20(Ig)%20are,)%2C%20and%2022%20(22q11.

ig_loci <- list(
  IGH = c(chr="14", start=105586437, end=106879844),
  IGK = c(chr="2",  start=88500000,  end=89300000),
  IGL = c(chr="22", start=22400000,  end=23700000)
)

is_in_ig <- function(chr, pos) {
  (chr == ig_loci$IGH["chr"] & pos >= as.numeric(ig_loci$IGH["start"]) & pos <= as.numeric(ig_loci$IGH["end"])) |
    (chr == ig_loci$IGK["chr"] & pos >= as.numeric(ig_loci$IGK["start"]) & pos <= as.numeric(ig_loci$IGK["end"])) |
    (chr == ig_loci$IGL["chr"] & pos >= as.numeric(ig_loci$IGL["start"]) & pos <= as.numeric(ig_loci$IGL["end"]))
}

muts_phase_filtered <- muts_phase %>% filter(!mapply(is_in_ig, chr, pos))

# Fit dndscv per group (overall ω)

groups <- c("WT","BL Pre-expansion","BL-Intermediate","BL-Private")

fit_by_group <- lapply(groups, function(g) {
  d <- filter(muts_phase_filtered, sampleID == g)
  if (nrow(d) == 0) return(NULL)
  dndscv(d %>% dplyr::select(sampleID, chr, pos, ref, mut),
         refdb = "hg38", outp = 1, max_coding_muts_per_sample = Inf)
})
names(fit_by_group) <- groups

dn_ds_overall <- do.call(rbind, Map(extract_overall, fit_by_group, names(fit_by_group))) %>%
  mutate(Group = factor(Group, levels = groups))

# Plot: dot + 95% CI

pal4 <- c("WT"="#D2BD96","BL Pre-expansion"="#1D3557",
          "BL-Intermediate"="#0A9086","BL-Private"="#A62639")

ggplot(dn_ds_overall, aes(Group, omega_all, color = Group)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, size = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = pal4) +
  scale_y_continuous(limits = c(0, 4), breaks = 0:4) +
  labs(title = "Overall dN/dS across coding classes",
       x = NULL, y = "dN/dS (ω)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")


# Order groups top-to-bottom

order_y <- c("BL-Private","BL-Intermediate", "BL Pre-expansion", "WT")

df_forest <- dn_ds_overall %>%
  mutate(Group = factor(Group, levels = order_y)) %>%
  arrange(Group)

# Plot as forest plot

p_forest <- ggplot(df_forest, aes(x = omega_all, y = Group)) +
  # horizontal CI
  geom_segment(aes(x = lo, xend = hi, yend = Group), linewidth = 0.9) +
  # point estimate (diamond)
  geom_point(shape = 18, size = 10) +
  # neutrality line
  geom_vline(xintercept = 1, linetype = "dashed") +
  # optional x-range (adjust to your data)
  scale_x_continuous(limits = c(0.5, 6.0), breaks = seq(0.5, 6.0, by = 0.5)) +
  labs(x = "dN/dS ratio (ω)", y = NULL, title = "Overall dN/dS per group") +
  theme_classic(base_size = 12) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/dNdS/dnds_forest_plot.pdf", plot = p_forest, width = 7, height = 6)

# add bulk WGS

vcf_paths <- c(
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P3G6/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P3G6_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PRN4/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PRN4_bulk.vep.SMuRF.filtered.sorted.vcf.gz",  
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P856/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P856_D.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/P856/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/P856_R.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PIA9/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PIA9_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PVA9/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PVA9_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/ASAP_FROM_CLOUD/PJBU/vcf_batches/batch_bulk/vcf/germline/somatic_filtering/SMuRF/PJBU_bulk.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID104AAO/SMuRF/PMCID104AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID132AAL/SMuRF/PMCID132AAL.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID137AAO/SMuRF/PMCID137AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID163AAN/SMuRF/PMCID163AAN.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID321AAO/SMuRF/PMCID321AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID340AAO/SMuRF/PMCID340AAO.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID491AAS/SMuRF/PMCID491AAS.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID509AAT/SMuRF/PMCID509AAT.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID540AAN/SMuRF/PMCID540AAN.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID610AAS/SMuRF/PMCID610AAS.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID690AAT/SMuRF/PMCID690AAT.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID821AAL/SMuRF/PMCID821AAL.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID867AAT/SMuRF/PMCID867AAT.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID967AAP/SMuRF/PMCID967AAP.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID211AAO_2/SMuRF/PMCID211AAO_2.vep.SMuRF.filtered.sorted.vcf.gz",
  "~/hpc/pmc_vanboxtel/projects/Burkitt/1_Input/Diagnostic_samples/PMCID458AAQ/SMuRF/PMCID458AAQ.vep.SMuRF.filtered.sorted.vcf.gz"
)

vcf_paths <- vcf_paths[!grepl("690AAT|458AAQ", vcf_paths)] # no myc translocation found in this sample

# --- Build Bulk WGS mutations from the 21 VCFs -------------------------------

stopifnot(length(vcf_paths) == 21)

# Extract sample names from VCF paths
samples_bulk <- basename(vcf_paths) |>
  gsub("\\.vep\\.SMuRF\\.filtered\\.sorted\\.vcf\\.gz$|\\.vcf\\.gz$", "", x = _)

# Function to convert one VCF to mutation table
vcf_to_muts_bulk <- function(vcf_path, sample_id) {
  v <- VariantAnnotation::readVcf(vcf_path, "hg38")
  rr <- rowRanges(v)
  alts <- unlist(VariantAnnotation::alt(v))
  
  data.frame(
    sampleID = sample_id,
    chr = as.character(GenomeInfoDb::seqnames(rr)),
    pos = BiocGenerics::start(rr),
    ref = as.character(VariantAnnotation::ref(v)),
    mut = as.character(alts),
    stringsAsFactors = FALSE
  )
}

# Apply to all bulk samples
muts_bulk <- dplyr::bind_rows(purrr::map2(vcf_paths, samples_bulk, vcf_to_muts_bulk)) %>%
  mutate(chr = as.character(chr)) %>%
  filter(chr %in% as.character(1:22)) %>%   # keep autosomes only
  filter(!mapply(is_in_ig, chr, pos))        # remove IGH, IGK, IGL regions

fit_bulk <- dndscv(muts_bulk, refdb = "hg38",
                   outp = 1, max_coding_muts_per_sample = Inf)

bulk_overall <- extract_overall(fit_bulk, "Bulk WGS")

dn_ds_overall2 <- dplyr::bind_rows(dn_ds_overall, bulk_overall) %>%
  dplyr::mutate(Group = factor(Group, levels = c("WT","BL Pre-expansion","BL-Intermediate","BL-Private","Bulk WGS")))

pal5 <- c("WT"="#D2BD96","BL Pre-expansion"="#1D3557",
          "BL-Intermediate"="#0A9086","BL-Private"="#A62639",
          "Bulk WGS"="#6C757D")

ggplot(dn_ds_overall2, aes(Group, omega_all, color = Group)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, size = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = pal5) +
  scale_y_continuous(limits = c(0, 4), breaks = 0:4) +
  labs(title = "Overall dN/dS across coding classes (incl. Bulk WGS)",
       x = NULL, y = "dN/dS (ω)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

order_y2 <- c("BL-Private","BL-Intermediate","BL Pre-expansion","Bulk WGS", "WT")

df_forest2 <- dn_ds_overall2 %>%
  dplyr::mutate(Group = factor(Group, levels = order_y2)) %>%
  dplyr::arrange(Group)

p_forest2 <- ggplot(df_forest2, aes(x = omega_all, y = Group)) +
  geom_segment(aes(x = lo, xend = hi, yend = Group), linewidth = 0.9) +
  geom_point(shape = 18, size = 10) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(limits = c(0.5, 6.0), breaks = seq(0.5, 6.0, by = 0.5)) +
  labs(x = "dN/dS ratio (ω)", y = NULL, title = "Overall dN/dS per group (incl. Bulk WGS)") +
  theme_classic(base_size = 12) +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

p_forest2

ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/dNdS/dnds_with_bulk_forest_plot.pdf", plot = p_forest2, width = 7, height = 6)

# --- Append Bulk WGS to the existing per-group results -----------------------

dn_ds_overall2 <- dplyr::bind_rows(dn_ds_overall, bulk_overall) %>%
  dplyr::mutate(Group = factor(Group, levels = c("WT","BL Pre-expansion","BL-Intermediate","BL-Private","Bulk WGS")))

# # Merge BL Pre-expansion (Trunk) + BL-Intermediate into a single group
# 
# muts_phase_merged <- muts_phase_filtered %>%
#   dplyr::mutate(sampleID = dplyr::case_when(
#     sampleID %in% c("BL Pre-expansion", "BL-Intermediate") ~ "BL Trunk+Intermediate",
#     TRUE ~ sampleID
#   ))
# 
# # Groups to fit
# 
# groups <- c("WT", "BL Trunk+Intermediate", "BL-Private")
# 
# # Fit dndscv per merged group
# 
# fit_by_group <- lapply(groups, function(g) {
#   d <- dplyr::filter(muts_phase_merged, sampleID == g)
#   if (nrow(d) == 0) return(NULL)
#   dndscv(d %>% dplyr::select(sampleID, chr, pos, ref, mut),
#          refdb = "hg38", outp = 1, max_coding_muts_per_sample = Inf)
# })
# names(fit_by_group) <- groups
# 
# # Extract overall dN/dS (ω)
# 
# dn_ds_overall <- do.call(
#   rbind,
#   Map(extract_overall, fit_by_group, names(fit_by_group))
# ) %>% dplyr::mutate(Group = factor(Group, levels = groups))
# 
# # --- Dot + 95% CI plot (3 groups) ---
# pal3 <- c("WT" = "#D2BD96",
#           "BL Trunk+Intermediate" = "#1D3557",
#           "BL-Private" = "#A62639")
# 
# ggplot(dn_ds_overall, aes(Group, omega_all, color = Group)) +
#   geom_point(size = 4) +
#   geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, size = 0.6) +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
#   scale_color_manual(values = pal3) +
#   scale_y_continuous(limits = c(0, 4), breaks = 0:4) +
#   labs(title = "Overall dN/dS across coding classes",
#        x = NULL, y = "dN/dS (ω)") +
#   theme_classic(base_size = 12) +
#   theme(legend.position = "none")
# 
# # Forest plot (3 groups) 
# 
# order_y <- c("BL-Private", "BL Trunk+Intermediate", "WT")
# df_forest <- dn_ds_overall %>%
#   dplyr::mutate(Group = factor(Group, levels = order_y)) %>%
#   dplyr::arrange(Group)
# 
# ggplot(df_forest, aes(x = omega_all, y = Group)) +
#   geom_segment(aes(x = lo, xend = hi, yend = Group), linewidth = 0.9) +
#   geom_point(shape = 18, size = 10) +
#   geom_vline(xintercept = 1, linetype = "dashed") +
#   scale_x_continuous(limits = c(0.5, 3.2), breaks = seq(0.5, 3.0, by = 0.5)) +
#   labs(x = "dN/dS ratio (ω)", y = NULL, title = "Overall dN/dS per group") +
#   theme_classic(base_size = 12) +
#   theme(axis.line.y = element_blank(),
#         axis.ticks.y = element_blank())
# 
# 
# 
# # Per-gene dN/dS
# 
# # POOLED per-gene dN/dS (uses already-filtered muts_phase_filtered)
# fit_all <- dndscv(
#   muts_phase_filtered %>% dplyr::select(sampleID, chr, pos, ref, mut),
#   refdb = "hg38", outp = 1, max_coding_muts_per_sample = Inf
# )
# 
# genes_all <- as_tibble(fit_all$genemuts)
# stopifnot(all(c("gene_name","n_syn","n_mis","n_non","n_spl","exp_syn","exp_mis","exp_non","exp_spl") %in% names(genes_all)))
# 
# # Robust Poisson test helper (obs vs expected under neutrality)
# 
# rate_test <- function(obs, exp) {
#   if (is.na(obs) || is.na(exp) || exp <= 0) {
#     return(c(omega = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_))
#   }
#   # poisson.test expects integer counts for x and a positive exposure T
#   pt <- poisson.test(x = as.integer(round(obs)), T = exp, r = 1)
#   c(
#     omega = obs / exp,
#     lo    = unname(pt$conf.int[1]),
#     hi    = unname(pt$conf.int[2]),
#     p     = unname(pt$p.value)
#   )
# }
# 
# calc_one <- function(n_mis, exp_mis, n_non, exp_non, n_spl, exp_spl) {
#   mis  <- rate_test(n_mis, exp_mis)
#   trnc <- rate_test(n_non, exp_non)
#   spl  <- rate_test(n_spl, exp_spl)
#   n_all   <- sum(n_mis, n_non, n_spl, na.rm = TRUE)
#   exp_all <- sum(exp_mis, exp_non, exp_spl, na.rm = TRUE)
#   wall <- rate_test(n_all, exp_all)
#   
#   c(
#     wmis      = mis["omega"],   wmis_lo   = mis["lo"],   wmis_hi   = mis["hi"],   p_mis   = mis["p"],
#     wtrunc    = trnc["omega"],  wtrunc_lo = trnc["lo"],  wtrunc_hi = trnc["hi"],  p_trunc = trnc["p"],
#     wspl      = spl["omega"],   wspl_lo   = spl["lo"],   wspl_hi   = spl["hi"],   p_spl   = spl["p"],
#     wall      = wall["omega"],  wall_lo   = wall["lo"],  wall_hi   = wall["hi"],  p_wall  = wall["p"]
#   )
# }
# 
# # Vectorised compute across all genes
# 
# calc <- genes_all %>%
#   transmute(n_mis, exp_mis, n_non, exp_non, n_spl, exp_spl) %>%
#   pmap_dfr(function(n_mis, exp_mis, n_non, exp_non, n_spl, exp_spl) {
#     out <- calc_one(n_mis, exp_mis, n_non, exp_non, n_spl, exp_spl)  # named numeric vector
#     as_tibble_row(as.list(out))  # coerce to 1-row tibble
#   })
# 
# # Generic name cleaner for your current headers
# 
# clean_names <- function(nms) {
#   nms <- sub("\\.omega$", "", nms)          # wmis.omega -> wmis
#   nms <- sub("_lo\\.lo$", "_lo", nms)       # wmis_lo.lo -> wmis_lo
#   nms <- sub("_hi\\.hi$", "_hi", nms)       # wmis_hi.hi -> wmis_hi
#   nms <- sub("^p_(.+)\\.p$", "p_\\1", nms)  # p_mis.p -> p_mis
#   nms
# }
# 
# calc_clean <- calc
# names(calc_clean) <- clean_names(names(calc_clean))
# 
# genes_omega <- dplyr::bind_cols(
#   genes_all %>% dplyr::select(gene = gene_name, n_syn:n_spl, exp_syn:exp_spl),
#   calc_clean
# ) %>%
#   dplyr::mutate(
#     n_nonsyn = n_mis + n_non + n_spl,
#     q_mis    = p.adjust(p_mis,   method = "BH"),
#     q_trunc  = p.adjust(p_trunc, method = "BH"),
#     q_spl    = p.adjust(p_spl,   method = "BH"),
#     q_wall   = p.adjust(p_wall,  method = "BH")
#   )
# 
# genes_omega_nonzero <- genes_omega %>%
#   filter(n_nonsyn > 0)
# 
# # Simple volcano (overall “wall”)
# 
# volc <- genes_omega_nonzero %>%
#   filter(!is.na(wall), !is.na(q_wall)) %>%
#   mutate(log10q = -log10(pmax(q_wall, 1e-300)),
#          sig = q_wall < 0.1)
# 
# toplab <- volc %>% arrange(q_wall) %>% slice_head(n = 20)
# 
# label_df <- volc %>%
#   filter(log10q > 0.5) %>%          # keep only genes above the cutoff
#   arrange(q_wall) %>%               # optional: rank by significance
#   slice_head(n = 20)                # label top 20 among them
# 
# p <- ggplot(volc, aes(x = wall, y = log10q)) +
#   geom_point(aes(alpha = sig)) +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
#   geom_hline(yintercept = -log10(0.10), linetype = "dotted", color = "grey50") +
#   ggrepel::geom_text_repel(
#     data = label_df,
#     aes(label = gene),
#     size = 3,
#     max.overlaps = 100
#   ) +
#   scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4)) +
#   labs(title = "Pooled per-gene dN/dS (overall)",
#        x = "ω (overall)",
#        y = expression(-log[10](q[FDR]))) +
#   theme_classic(base_size = 12) +
#   theme(legend.position = "none")
# 
# # save as PDF (adjust width/height in inches)
# ggsave("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/CellPhy/Figures/dNdS/dnds_volcano.pdf", plot = p, width = 7, height = 6)
# 
