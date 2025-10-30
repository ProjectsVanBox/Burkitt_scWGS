# ============================
# BL signatures: trunk vs post-expansion (faceted, with p-values)
# ============================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(lme4)
library(emmeans)
library(glmmTMB)

# ---- Palette (use yours if already defined) ----
if (!exists("pal")) {
  pal <- c("#D2BD96", "#0A9086", "#B3B3B3", "#A62639", "#1D3557")
  names(pal) <- c("SBS1","SBS9","SBS17b","SBS18","SBSblood")
}

# ---- 1) Read files (tab or comma) ----
dir_path <- "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Data"

files <- list.files(dir_path, pattern = "^sbs_contr_per_group_.*\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)

get_pid <- function(f) sub("^sbs_contr_per_group_(.*)\\.csv$", "\\1", basename(f))

patient_mats <- setNames(lapply(files, function(f) {
  read.csv(f, header = TRUE, row.names = 1, check.names = FALSE)
}) |> lapply(as.matrix),
sapply(files, get_pid))



# signatures (rows we keep/order)
sig_keep  <- c("SBS1","SBS17b","SBS18","SBS9","SBSblood")
# groups (columns we compare/order)
keep_cols <- c("BL_Trunk","BL_post_expansion")

# 1) Prep per-patient matrices: drop PTA, keep rows/cols we need
patient_mats2 <- lapply(names(patient_mats), function(pid) {
  mat <- patient_mats[[pid]]
  if ("PTA_v1" %in% rownames(mat)) {
    mat <- mat[setdiff(rownames(mat), "PTA_v1"), , drop = FALSE]
  }
  miss_cols <- setdiff(keep_cols, colnames(mat))
  if (length(miss_cols)) stop(sprintf("Patient %s missing: %s", pid, paste(miss_cols, collapse=", ")))
  miss_sig  <- setdiff(sig_keep, rownames(mat))
  if (length(miss_sig))  stop(sprintf("Patient %s missing signature rows: %s", pid, paste(miss_sig, collapse=", ")))
  mat[sig_keep, keep_cols, drop = FALSE]
})
names(patient_mats2) <- names(patient_mats)

# 2) Long format + relative contribution within patient & group
df_long <- bind_rows(lapply(names(patient_mats2), function(pid) {
  m <- patient_mats2[[pid]]
  data.frame(Patient = pid,
             Signature = rownames(m),
             BL_Trunk = m[, "BL_Trunk"],
             BL_post_expansion = m[, "BL_post_expansion"],
             check.names = FALSE)
})) |>
  pivot_longer(cols = c(BL_Trunk, BL_post_expansion),
               names_to = "Group", values_to = "Absolute") |>
  group_by(Patient, Group) |>
  mutate(Relative = Absolute / sum(Absolute)) |>
  ungroup() |>
  mutate(
    Group = factor(Group, levels = keep_cols),
    Signature = factor(Signature, levels = sig_keep)
  )


# 3) Mixed-effects EMMs per signature (beta GLMM) + CI + p-values (two-sided, unadjusted)
eps <- 1e-6
df_beta <- df_long |>
  mutate(Rel_beta = (Relative * (1 - 2*eps)) + eps)  # push into (0,1) for beta

emm_list <- list(); p_list <- list()

for (sig in levels(df_beta$Signature)) {
  d <- df_beta %>% filter(Signature == sig)
  
  fit <- glmmTMB(Rel_beta ~ Group + (1|Patient),
                 data = d, family = beta_family(link = "logit"))
  
  # EMMs on RESPONSE scale with 95% CI
  emms <- emmeans(fit, ~ Group)
  # Symmetric (delta-method) CIs on the response scale
  link_sum <- as.data.frame(summary(emms, type = "link", infer = c(TRUE, TRUE)))
  # link_sum has columns: Group, emmean (η), SE, lower.CL, upper.CL on the logit scale
  
  eta   <- link_sum$emmean
  seEta <- link_sum$SE
  p     <- plogis(eta)            # back-transform mean to 0–1
  gprime <- p * (1 - p)           # derivative of invlogit at η
  
  seP   <- seEta * gprime         # delta-method SE on probability scale
  lower_sym <- p - 1.96 * seP     # symmetric 95% CI around p
  upper_sym <- p + 1.96 * seP
  
  emm_df <- data.frame(
    Group    = link_sum$Group,
    emmean   = p,
    lower    = pmax(0, lower_sym),  # clamp to [0,1]
    upper    = pmin(1, upper_sym),
    Signature = sig
  )
  
  # Two-sided contrast: (Post − Trunk) on response scale; unadjusted p
  diff_df <- as.data.frame(contrast(emms, list("Post - Trunk" = c(-1, 1)), type = "response"))
  p_list[[sig]] <- diff_df %>%
    transmute(Signature = sig, p_value = p.value)
  
  emm_list[[sig]] <- emm_df
}

emm_all <- dplyr::bind_rows(emm_list) %>%
  dplyr::mutate(
    Signature = factor(Signature, levels = sig_keep),
    Group     = factor(as.character(Group), levels = keep_cols),
    # keep everything within [0,1] just in case of tiny numeric drift
    emmean = pmin(pmax(emmean, 0), 1),
    lower  = pmin(pmax(lower , 0), 1),
    upper  = pmin(pmax(upper , 0), 1)
  )

pvals <- dplyr::bind_rows(p_list) %>%
  dplyr::mutate(
    Signature = factor(Signature, levels = sig_keep),
    p_label   = paste0("p = ", formatC(p_value, format = "e", digits = 2))
  )

# Annotation positions (still 0–1)
annot <- emm_all %>%
  dplyr::group_by(Signature) %>%
  dplyr::summarise(y_pos = pmin(1, max(upper, na.rm = TRUE) * 1.08), .groups = "drop") %>%
  dplyr::left_join(pvals, by = "Signature") %>%
  dplyr::mutate(x_pos = 1.5)

# ---- Plot on 0–1 axis ----
p <- ggplot(emm_all, aes(x = Group, y = emmean, color = Signature)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.18) +
  facet_wrap(~ Signature, nrow = 1) +
  scale_color_manual(values = pal, drop = FALSE, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    x = NULL,
    y = "Relative contribution",
    title = "Cohort-wide estimated proportions: Pre-expansion vs Post-expansion",
    subtitle = "The dots represent the mixed effect meta-analysis estimated proportion and the bars 95% confidence intervals. The P values are for a two-sided test."
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  ) +
  geom_text(data = annot, aes(x = x_pos, y = y_pos, label = p_label),
            inherit.aes = FALSE, size = 3.2)

print(p)

# Save as PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/MutationalPatterns/Figures/Pre_vs_Post_expansion_sig_comparison.pdf",
  plot = p,
  width = 8,
  height = 6
)
