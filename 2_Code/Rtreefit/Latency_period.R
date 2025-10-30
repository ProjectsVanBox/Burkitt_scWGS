library(readxl)
library(ggplot2)

df <- read_excel("~/surfdrive/Shared/pmc_vanboxtel/projects/Burkitt_github/3_Output/Rtreefit/Data/latency_period.xlsx")

df_filtered <- df %>% filter(Patient != "P856_R")

ggplot(df_filtered, aes(x = sampling_age, y = `latency period`, label = Patient)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(vjust = -0.5, size = 3) +
  geom_line(group = 1, color = "grey50", linetype = "dashed") +
  labs(
    x = "Sampling age (years)",
    y = "Latency period (years)",
    title = "Latency period vs. Sampling age"
  ) +
  theme_minimal(base_size = 14)

ggplot(df %>% filter(Patient != "P856_R"),
       aes(x = sampling_age, y = `latency period`, label = Patient)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(vjust = -0.5, size = 3) +
  geom_smooth(method = "loess", se = TRUE, color = "darkgreen") +
  labs(
    x = "Sampling age (years)",
    y = "Latency period (years)",
    title = "Latency period vs. Sampling age (loess smooth)"
  ) +
  theme_minimal(base_size = 14)
