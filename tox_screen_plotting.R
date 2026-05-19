#tox screen

# ================================
# Toxicity screen Try 3 analysis
# Dox and TMP temporal plots
# ================================

library(tidyverse)
library(readr)
source("nathan_functs.R")
output_dir <- "final_images"

# ----------------
# 1. Read data
# ----------------

tox <- read_csv(
  "Luminescence tox screen 2 and 3 - Try 3.csv",
  col_types = cols(.default = col_character())
)

# Rename columns manually because the luminescence column is unnamed
colnames(tox) <- c("time_hr", "dox_ug_ml", "tmp_uM", "replicate", "signal")

# ----------------
# 2. Clean data
# ----------------

# Molecular weight for doxycycline hyclate
doxy_mw <- 512.94

tox_clean <- tox %>%
  mutate(
    time_hr = as.numeric(time_hr),
    dox_ug_ml = na_if(dox_ug_ml, "Blank"),
    tmp_uM = na_if(tmp_uM, "Blank"),
    dox_ug_ml = as.numeric(dox_ug_ml),
    tmp_uM = as.numeric(tmp_uM),
    replicate = as.numeric(replicate),
    signal = parse_number(signal, locale = locale(grouping_mark = ",")),
    
    # Convert Dox from ug/mL to uM
    # ug/mL = mg/L, so uM = ug/mL * 1000 / molecular weight
    dox_uM = dox_ug_ml * 1000 / doxy_mw
  )
# ----------------
# Blank correction
# ----------------

blank_by_time <- tox_clean %>%
  filter(is.na(dox_ug_ml), is.na(tmp_uM)) %>%
  group_by(time_hr) %>%
  summarise(
    blank_signal = mean(signal, na.rm = TRUE),
    .groups = "drop"
  )

tox_norm <- tox_clean %>%
  left_join(blank_by_time, by = "time_hr") %>%
  mutate(
    blank_corrected_signal = signal - blank_signal
  )

# ----------------
# Dox growth-corrected relative viability
# Each condition starts at 100% at 0 h
# Dox 0 stays around 100% over time
# TMP = 0
# ----------------

dox_data <- tox_norm %>%
  filter(
    tmp_uM == 0,
    !is.na(dox_ug_ml),
    time_hr <= 72
  )

# Mean signal at 0 h for each Dox concentration
dox_baseline <- dox_data %>%
  filter(time_hr == 0) %>%
  group_by(dox_uM) %>%
  summarise(
    baseline_signal = mean(blank_corrected_signal, na.rm = TRUE),
    .groups = "drop"
  )

# Normalize each concentration to its own 0 h
dox_growth <- dox_data %>%
  left_join(dox_baseline, by = "dox_uM") %>%
  mutate(
    growth_relative_to_0h = blank_corrected_signal / baseline_signal
  )

# Get the untreated Dox 0 growth curve
dox_control_growth <- dox_growth %>%
  filter(dox_ug_ml == 0) %>%
  group_by(time_hr) %>%
  summarise(
    control_growth = mean(growth_relative_to_0h, na.rm = TRUE),
    .groups = "drop"
  )

# Normalize each condition's growth to untreated growth at the same time point
dox_summary <- dox_growth %>%
  left_join(dox_control_growth, by = "time_hr") %>%
  mutate(
    relative_viability = 100 * growth_relative_to_0h / control_growth,
    relative_viability = pmax(relative_viability, 0)  # floor negatives only
  ) %>%
  group_by(time_hr, dox_ug_ml, dox_uM) %>%
  summarise(
    mean_viability = mean(relative_viability, na.rm = TRUE),
    sd_viability = sd(relative_viability, na.rm = TRUE),
    n = sum(!is.na(relative_viability)),
    se_viability = sd_viability / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    # Plot-only capped values
    mean_viability_plot = pmin(mean_viability, 100),
    ymin_plot = pmax(mean_viability - se_viability, 0),
    ymax_plot = pmin(mean_viability + se_viability, 100),
    
    dox_label = case_when(
      dox_ug_ml == 0 ~ "0 µM",
      dox_ug_ml == 1 ~ "2 µM",
      dox_ug_ml == 10 ~ "20 µM",
      dox_ug_ml == 50 ~ "100 µM",
      dox_ug_ml == 100 ~ "200 µM",
      TRUE ~ paste0(round(dox_uM, 0), " µM")
    ),
    dox_label = factor(
      dox_label,
      levels = c("0 µM", "2 µM", "20 µM", "100 µM", "200 µM")
    )
  )

dox_plot <- ggplot(
  dox_summary,
  aes(
    x = time_hr,
    y = mean_viability,
    color = dox_label,
    group = dox_label
  )
) +
  geom_hline(
    yintercept = 100,
    linetype = "dotted",
    linewidth = 0.6,
    color = "gray40"
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(
      ymin = mean_viability - se_viability,
      ymax = mean_viability + se_viability
    ),
    width = 2,
    linewidth = 0.4
  ) +
  labs(
    x = "Time (hr)",
    y = "Relative cell viability",
    color = "Dox concentration"
  ) +
  scale_x_continuous(breaks = c(0, 24, 48, 72)) +
  coord_cartesian(ylim = c(0, 120)) +
  scale_color_manual(
    values = c(
      "0 µM" = "#b4d9cc",
      "2 µM" = "#89c0b6",
      "20 µM" = "#63a6a0",
      "100 µM" = "#448c8a",
      "200 µM" = "#287274"
    )
  ) +
  theme_nathan()

print(dox_plot)

ggsave(file.path(output_dir,"Tox_Screen_Dox_Try3_Internal_Normalized.png"), dox_plot, width = 4, height = 2.5, dpi = 300)

ggsave(file.path(output_dir,"Tox_Screen_Dox_Try3_Internal_Normalized.pdf"), dox_plot, width = 4, height = 2.5)


# ----------------
# TMP growth-corrected relative viability
# Each condition starts at 100% at 0 h
# TMP 0 stays around 100% over time
# Dox = 0
# ----------------

tmp_data <- tox_norm %>%
  filter(
    dox_ug_ml == 0,
    !is.na(tmp_uM),
    time_hr <= 72
  )

# Mean signal at 0 h for each TMP concentration
tmp_baseline <- tmp_data %>%
  filter(time_hr == 0) %>%
  group_by(tmp_uM) %>%
  summarise(
    baseline_signal = mean(blank_corrected_signal, na.rm = TRUE),
    .groups = "drop"
  )

# Normalize each concentration to its own 0 h
tmp_growth <- tmp_data %>%
  left_join(tmp_baseline, by = "tmp_uM") %>%
  mutate(
    growth_relative_to_0h = blank_corrected_signal / baseline_signal
  )

# Get the untreated TMP 0 growth curve
tmp_control_growth <- tmp_growth %>%
  filter(tmp_uM == 0) %>%
  group_by(time_hr) %>%
  summarise(
    control_growth = mean(growth_relative_to_0h, na.rm = TRUE),
    .groups = "drop"
  )

# Normalize each condition's growth to untreated growth at the same time point
tmp_summary <- tmp_growth %>%
  left_join(tmp_control_growth, by = "time_hr") %>%
  mutate(
    relative_viability = 100 * growth_relative_to_0h / control_growth,
    relative_viability = pmax(relative_viability, 0)   # floor negatives only
  ) %>%
  group_by(time_hr, tmp_uM) %>%
  summarise(
    mean_viability = mean(relative_viability, na.rm = TRUE),
    sd_viability = sd(relative_viability, na.rm = TRUE),
    n = sum(!is.na(relative_viability)),
    se_viability = sd_viability / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    # Optional plot-only capped values
    mean_viability_plot = pmin(mean_viability, 100),
    ymin_plot = pmax(mean_viability - se_viability, 0),
    ymax_plot = pmin(mean_viability + se_viability, 100),
    
    tmp_label = factor(
      paste0(tmp_uM, " µM"),
      levels = paste0(c(0, 0.1, 1, 10, 100), " µM")
    )
  )

tmp_plot <- ggplot(
  tmp_summary,
  aes(
    x = time_hr,
    y = mean_viability_plot,
    color = tmp_label,
    group = tmp_label
  )
) +
  geom_hline(
    yintercept = 100,
    linetype = "dotted",
    linewidth = 0.6,
    color = "gray40"
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(
      ymin = ymin_plot,
      ymax = ymax_plot
    ),
    width = 2,
    linewidth = 0.4
  ) +
  labs(
    x = "Time (hr)",
    y = "Relative cell viability",
    color = "TMP concentration"
  ) +
  scale_x_continuous(breaks = c(0, 24, 48, 72)) +
  coord_cartesian(ylim = c(0, 120)) +
  scale_color_manual(
    values = c(
      "0 µM" = "#b4d9cc",
      "0.1 µM" = "#89c0b6",
      "1 µM" = "#63a6a0",
      "10 µM" = "#448c8a",
      "100 µM" = "#287274"
    )
  ) +
  theme_nathan()

print(tmp_plot)

ggsave(file.path(output_dir,"Tox_Screen_TMP_Try3_Internal_Normalized.png"), tmp_plot, width = 4, height = 2.5, dpi = 300)
ggsave(file.path(output_dir,"Tox_Screen_TMP_Try3_Internal_Normalized.pdf"), tmp_plot, width = 4, height = 2.5)
