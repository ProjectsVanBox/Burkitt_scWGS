################################################################################
# Manuscript: Clonal Evolution of Paediatric Burkitt Lymphoma Through Time and Space
# Description: latency correlations
# Author: Alexander Steemers
################################################################################

# Load libraries

library(readxl)
library(dplyr)
library(ggplot2)
source('~/hpc/pmc_vanboxtel/projects/CHemALL/2_Code/theme_CHemALL.R')

# Load data

df <- read_excel("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/1_Input/latency_period.xlsx")

# Exponential model (log-linear)

df <- df %>% mutate(log_mut_rate = log(Mut_rate_tumour_mean))
exp_lm <- lm(log_mut_rate ~ Latency_period, data = df)

# Stats for annotation

r2_exp  <- summary(exp_lm)$r.squared
p_val   <- summary(exp_lm)$coefficients["Latency_period", "Pr(>|t|)"]
p_str   <- ifelse(p_val < 1e-3, "< 0.001", sprintf("%.3f", p_val))

# Smooth predictions (for curve + 95% CI), over plotting range

x_min <- 1
x_max <- 10
newdata <- data.frame(Latency_period = seq(x_min, x_max, length.out = 200))
pred    <- predict(exp_lm, newdata, interval = "confidence")

curve_df <- newdata %>%
  mutate(
    fit = exp(pred[, "fit"]),
    lwr = exp(pred[, "lwr"]),
    upr = exp(pred[, "upr"])
  )

# Axes

y_min <- 0
y_max <- 350

# Annotation placement (inside the panel)

ann_x <- x_max - 2
ann_y <- y_max * 0.80

# Plot

p1 <- ggplot(df, aes(x = Latency_period, y = Mut_rate_tumour_mean)) +
  geom_point(size = 3, color = "black", na.rm = TRUE) +
  
  # 95% confidence ribbon of exponential model
  geom_ribbon(
    data = curve_df,
    aes(x = Latency_period, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +
  
  # Exponential regression curve
  geom_line(
    data = curve_df,
    aes(x = Latency_period, y = fit),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1
  ) +
  
  annotate(
    "text",
    x = ann_x, y = ann_y, hjust = 0, size = 5,
    label = paste0("R² = ", sprintf("%.3f", r2_exp), "\n",
                   "p = ", p_str)
  ) +
  labs(
    x = "Latency period (years)",
    y = "Mean tumour mutation rate post-expansion (SBS/year)"
  ) +
  scale_x_continuous(limits = c(x_min, x_max)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = seq(y_min, y_max, 50)) +
  theme_CHemALL() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

print(p1)

# Save as PDF

ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/mut_rate_vs_latency.pdf",
  plot = p1,
  width = 8, height = 6
)

# Exponential model (log-linear): log(mut rate) ~ Age
df <- df %>% mutate(log_mut_rate = log(Mut_rate_tumour_mean))
exp_lm_age <- lm(log_mut_rate ~ Age, data = df)

# Stats for annotation
r2_exp_age  <- summary(exp_lm_age)$r.squared
p_age       <- summary(exp_lm_age)$coefficients["Age", "Pr(>|t|)"]
p_str_age   <- ifelse(p_age < 1e-3, "< 0.001", sprintf("%.3f", p_age))

# Prediction grid over Age range
x_min_age <- min(df$Age, na.rm = TRUE)
x_max_age <- max(df$Age, na.rm = TRUE)
newdata_age <- data.frame(Age = seq(x_min_age, x_max_age, length.out = 200))
pred_age    <- predict(exp_lm_age, newdata_age, interval = "confidence")

curve_df_age <- newdata_age %>%
  mutate(
    fit = exp(pred_age[, "fit"]),
    lwr = exp(pred_age[, "lwr"]),
    upr = exp(pred_age[, "upr"])
  )

# Axes (same y-scale as before)
y_min <- 0
y_max <- 350

# Annotation placement
ann_x_age <- x_max_age - 2
ann_y_age <- y_max * 0.80

# Plot
p_age <- ggplot(df, aes(x = Age, y = Mut_rate_tumour_mean)) +
  geom_point(size = 3, color = "black", na.rm = TRUE) +
  geom_ribbon(
    data = curve_df_age,
    aes(x = Age, ymin = lwr, ymax = upr),
    inherit.aes = FALSE, alpha = 0.20
  ) +
  geom_line(
    data = curve_df_age,
    aes(x = Age, y = fit),
    inherit.aes = FALSE, color = "black", linewidth = 1
  ) +
  annotate(
    "text",
    x = ann_x_age, y = ann_y_age, hjust = 0, size = 5,
    label = paste0("R² = ", sprintf("%.3f", r2_exp_age), "\n",
                   "p = ", p_str_age)
  ) +
  labs(
    x = "Age at sampling (years)",
    y = "Mean tumour mutation rate post-expansion (SBS/year)"
  ) +
  scale_x_continuous(limits = c(x_min_age, x_max_age)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = seq(y_min, y_max, 50)) +
  theme_CHemALL() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

print(p_age)

# Save PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/mut_rate_post_MRCA_sample_age.pdf",
  plot = p_age,
  width = 8, height = 6
)

# Slope + CI from the log-linear model
b      <- coef(exp_lm_age)["Age"]
ci     <- confint(exp_lm_age)["Age", ]
pct    <- (exp(b)  - 1) * 100
pct_ci <- (exp(ci) - 1) * 100

slope_txt <- paste0(
  "Slope (log-scale) β = ", sprintf("%.3f", b),
  "  [", sprintf("%.3f", ci[1]), ", ", sprintf("%.3f", ci[2]), "]",
  "\n≈ ", sprintf("%.1f", pct), "% / year  [",
  sprintf("%.1f", pct_ci[1]), ", ", sprintf("%.1f", pct_ci[2]), "]"
)

p_age2 <- p_age +
  annotate("text",
           x = ann_x_age, y = ann_y_age * 0.9, hjust = 0, size = 5,
           label = slope_txt
  )

print(p_age2)

# Linear model: mutation rate ~ Age
lin_lm_age <- lm(Mut_rate_tumour_mean ~ Age, data = df)

# Stats for annotation
r2_lin_age <- summary(lin_lm_age)$r.squared
p_lin_age  <- summary(lin_lm_age)$coefficients["Age", "Pr(>|t|)"]
p_str_lin_age <- ifelse(p_lin_age < 1e-3, "< 0.001", sprintf("%.3f", p_lin_age))

# Prediction grid over Age range
x_min_age <- min(df$Age, na.rm = TRUE)
x_max_age <- max(df$Age, na.rm = TRUE)
newdata_age_lin <- data.frame(Age = seq(x_min_age, x_max_age, length.out = 200))
pred_age_lin    <- predict(lin_lm_age, newdata_age_lin, interval = "confidence")
curve_df_age_lin <- cbind(newdata_age_lin, as.data.frame(pred_age_lin))

# Axes
y_min <- 0
y_max <- 350

# Annotation placement
ann_x_age <- x_max_age - 2
ann_y_age <- y_max * 0.80

# Plot (linear)
p_age_lin <- ggplot(df, aes(x = Age, y = Mut_rate_tumour_mean)) +
  geom_point(size = 3, color = "black", na.rm = TRUE) +
  
  # 95% confidence ribbon
  geom_ribbon(
    data = curve_df_age_lin,
    aes(x = Age, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +
  
  # Linear regression line
  geom_line(
    data = curve_df_age_lin,
    aes(x = Age, y = fit),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1
  ) +
  
  # Annotation
  annotate(
    "text",
    x = ann_x_age, y = ann_y_age, hjust = 0, size = 5,
    label = paste0("R² = ", sprintf("%.3f", r2_lin_age), "\n",
                   "p = ", p_str_lin_age)
  ) +
  
  labs(
    x = "Age at sampling (years)",
    y = "Mean tumour mutation rate post-expansion (SBS/year)"
  ) +
  scale_x_continuous(limits = c(x_min_age, x_max_age)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = seq(y_min, y_max, 50)) +
  theme_CHemALL() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

print(p_age_lin)

# Extract slope (Age coefficient)
coef(lin_lm_age)["Age"]
# or more explicitly
slope <- summary(lin_lm_age)$coefficients["Age", "Estimate"]
slope

# Linear model
lin_lm <- lm(Mut_rate_tumour_mean ~ Latency_period, data = df)

# Stats for annotation
r2_lin <- summary(lin_lm)$r.squared
p_lin  <- summary(lin_lm)$coefficients["Latency_period", "Pr(>|t|)"]
p_str_lin <- ifelse(p_lin < 1e-3, "< 0.001", sprintf("%.3f", p_lin))

# Predictions over plotting range
newdata_lin <- data.frame(Latency_period = seq(x_min, x_max, length.out = 200))
pred_lin    <- predict(lin_lm, newdata_lin, interval = "confidence")
curve_lin   <- cbind(newdata_lin, as.data.frame(pred_lin))

# Plot (linear)
p_lin <- ggplot(df, aes(x = Latency_period, y = Mut_rate_tumour_mean)) +
  geom_point(size = 3, color = "black", na.rm = TRUE) +
  
  # 95% confidence ribbon
  geom_ribbon(
    data = curve_lin,
    aes(x = Latency_period, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.20
  ) +
  
  # Linear regression line
  geom_line(
    data = curve_lin,
    aes(x = Latency_period, y = fit),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1
  ) +
  
  # Annotation
  annotate(
    "text",
    x = ann_x, y = ann_y, hjust = 0, size = 5,
    label = paste0("R² = ", sprintf("%.3f", r2_lin), "\n",
                   "p = ", p_str_lin)
  ) +
  
  labs(
    x = "Latency period (years)",
    y = "Mean tumour mutation rate post-expansion (SBS/year)"
  ) +
  scale_x_continuous(limits = c(x_min, x_max)) +
  scale_y_continuous(limits = c(y_min, y_max), breaks = seq(y_min, y_max, 50)) +
  theme_CHemALL() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

print(p_lin)

# Fit exponential model: log(Latency)

df <- df %>% mutate(log_latency = log(Latency_period))
exp_lm_age <- lm(log_latency ~ Age, data = df)

# Stats for annotation

r2_age  <- summary(exp_lm_age)$r.squared
p_age   <- summary(exp_lm_age)$coefficients["Age", "Pr(>|t|)"]
p_age_str <- ifelse(p_age < 1e-3, "< 0.001", sprintf("%.3f", p_age))

# Prediction range

x_min_age <- min(df$Age, na.rm = TRUE)
x_max_age <- max(df$Age, na.rm = TRUE)

new_age <- data.frame(Age = seq(x_min_age, x_max_age, length.out = 200))
pred_age <- predict(exp_lm_age, new_age, interval = "confidence")

curve_df_age <- new_age %>%
  mutate(
    fit = exp(pred_age[, "fit"]),
    lwr = exp(pred_age[, "lwr"]),
    upr = exp(pred_age[, "upr"])
  )

# Y-axis limits for latency

y_min_lat <- 0
y_max_lat <- max(df$Latency_period) * 2

# Annotation placement

ann_x_age <- x_max_age - 2
ann_y_age <- y_max_lat * 0.85

# Plot

p2 <- ggplot(df, aes(x = Age, y = Latency_period)) +
  geom_point(size = 3, color = "black", na.rm = TRUE) +
  
  geom_ribbon(
    data = curve_df_age,
    aes(x = Age, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = curve_df_age,
    aes(x = Age, y = fit),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1
  ) +
  annotate("text",
           x = ann_x_age, y = ann_y_age, hjust = 0, size = 5,
           label = paste0("R² = ", sprintf("%.3f", r2_age), "\n",
                          "p = ", p_age_str)
  ) +
  labs(
    x = "Age at sampling (years)",
    y = "Latency period (years)"
  ) +
  scale_x_continuous(limits = c(x_min_age, x_max_age)) +
  scale_y_continuous(limits = c(y_min_lat, y_max_lat)) +
  theme_CHemALL() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

print(p2)

# Save PDF
ggsave(
  filename = "~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Figures/latency_vs_age.pdf",
  plot = p2,
  width = 8, height = 6
)


# Linear model: Latency ~ Age
lin_lm_age <- lm(Latency_period ~ Age, data = df)

# Stats
r2_age_lin <- summary(lin_lm_age)$r.squared
p_age_lin  <- summary(lin_lm_age)$coefficients["Age", "Pr(>|t|)"]
p_age_lin_str <- ifelse(p_age_lin < 1e-3, "< 0.001", sprintf("%.3f", p_age_lin))

# Predictions over plotting range
new_age_lin <- data.frame(Age = seq(x_min_age, x_max_age, length.out = 200))
pred_age_lin <- predict(lin_lm_age, new_age_lin, interval = "confidence")
curve_df_age_lin <- cbind(new_age_lin, as.data.frame(pred_age_lin))

# Plot (linear)
p2_lin <- ggplot(df, aes(x = Age, y = Latency_period)) +
  geom_point(size = 3, color = "black", na.rm = TRUE) +
  geom_ribbon(
    data = curve_df_age_lin,
    aes(x = Age, ymin = lwr, ymax = upr),
    inherit.aes = FALSE, alpha = 0.2
  ) +
  geom_line(
    data = curve_df_age_lin,
    aes(x = Age, y = fit),
    inherit.aes = FALSE, color = "black", linewidth = 1
  ) +
  annotate(
    "text",
    x = ann_x_age, y = ann_y_age, hjust = 0, size = 5,
    label = paste0("R² = ", sprintf("%.3f", r2_age_lin), "\n",
                   "p = ", p_age_lin_str)
  ) +
  labs(
    x = "Age at sampling (years)",
    y = "Latency period (years)"
  ) +
  scale_x_continuous(limits = c(x_min_age, x_max_age)) +
  scale_y_continuous(limits = c(y_min_lat, y_max_lat)) +
  theme_CHemALL() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

print(p2_lin)

AIC(lin_lm, exp_lm) 

# AIC differences (ΔAIC) greater than 10 are considered decisive evidence in favour of the model with the lower AIC.
# Here, the ΔAIC ≈ 75, which is huge — a very strong indication that the exponential model fits your data far better.


df$Patient_Age <- paste(df$Patient, df$Age, sep = "_")


df$Patient_Age <- factor(df$Patient_Age, levels = df$Patient_Age)

p3 <- ggplot(df, aes(x = Patient_Age, y = Mean_post_MRCA_mutation_load)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme_CHemALL() +
  labs(
    x = "Patient Age",
    y = "Mean post-MRCA Mutation Load",
    title = "Mutation Load post-MRCA by Patient Age"
  )

print(p3)

p4 <- ggplot(df, aes(x = Latency_period, y = Mean_post_MRCA_mutation_load)) +
  geom_point(size = 3, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  theme_minimal() +
  theme_CHemALL() +
  labs(
    x = "Latency Period",
    y = "Mean post-MRCA Mutation Load",
    title = "Relationship Between Latency Period and Mutation Load"
  )

p4

p4 <- ggplot(df, aes(x = Age, y = Mean_post_MRCA_mutation_load)) +
  geom_point(size = 3, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, color = "darkred") +
  theme_minimal() +
  theme_CHemALL() +
  labs(
    x = "Latency Period",
    y = "Mean post-MRCA Mutation Load",
    title = "Relationship Between Latency Period and Mutation Load"
  )

p4

