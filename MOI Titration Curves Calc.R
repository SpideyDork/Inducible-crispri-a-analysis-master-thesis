library(tidyverse)

# 1. Input the New Data from your latest results
# vol_n1 corresponds to -1, vol_n2 to -2, etc.
df_raw <- tribble(
  ~condition,       ~target_gene, ~vol_10, ~vol_1, ~vol_n1, ~vol_n2, ~vol_n3, ~vol_n4, ~vol_n5, ~no_virus, ~no_puro,
  "CRISPRa - Avg",  "NTC",        668996,  243813, 178521,  126050,  177319,  184513,  173108,  142468,    893992,
  "CRISPRa - Avg",  "GFI1B",      884630,  628428, 234485,  91370,   153714,  191609,  194581,  139889,    795311,
  "CRISPRa - Avg",  "CD46",       837778,  439794, 204689,  108343,  213955,  170925,  160970,  176938,    919826,
  "CRISPRi - Avg",  "NTC",        734014,  534373, 489477,  475161,  481399,  498251,  466961,  404751,    922417,
  "CRISPRi - Avg",  "GFI1B",      723563,  757642, 549082,  482399,  446301,  457066,  427969,  367587,    795853,
  "CRISPRi - Avg",  "CD46",       928353,  639479, 500018,  425978,  472693,  466228,  448588,  444214,    801307
)

# 2. Volume mapping logic (converts labels to uL)
vol_map <- c(
  "vol_10" = 10, "vol_1" = 1, "vol_n1" = 0.1, "vol_n2" = 0.01, 
  "vol_n3" = 0.001, "vol_n4" = 0.0001, "vol_n5" = 0.00001
)

# 3. Process data and calculate MOI
tdata_moi <- df_raw %>%
  # Split 'CRISPRa - Avg' into 'CRISPRa' and 'Avg'
  separate(condition, into = c("system", "replicate"), sep = " - ") %>%
  pivot_longer(
    cols = starts_with("vol_"),
    names_to = "volume_label",
    values_to = "cell_count"
  ) %>%
  mutate(
    volume_ul = vol_map[volume_label],
    # Survival = (Count - Background) / (TotalPossible - Background)
    pct_survival = (cell_count - no_virus) / (no_puro - no_virus) * 100,
    pct_survival = pmax(0, pmin(100, pct_survival)),
    # MOI = -ln(1 - fraction_infected)
    moi = -log(1 - (pct_survival / 100) * 0.999) # 0.999 to avoid Infinity at 100%
  )

# 4. View results specifically for your 0.2 - 0.4 range
print("Conditions near the 0.2 - 0.4 MOI range:")
tdata_moi %>%
  filter(volume_ul == 1) %>% # Looking at the 1uL volume where most data sits
  select(system, target_gene, volume_ul, pct_survival, moi) %>%
  arrange(moi)

#Generate plots
axis_breaks <- c(10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001)

axis_labels <- c(
  "10",
  "1",
  "10^-1",
  "10^-2",
  "10^-3",
  "10^-4",
  "10^-5"
)



plot_moi <- function(sys_name) {
  target_colors <- if (sys_name == "CRISPRi") {
    c(
      "NTC" = "#f4a3a8",
      "GFI1B" = "#8b3058",
      "CD46" = "#cc607d"
    )
  } else {
    c(
      "NTC" = "darkseagreen2",
      "GFI1B" = "darkseagreen4",
      "CD46" = "darkseagreen3"
    )
  }
  #crispri
  #ffc6c4,#f4a3a8,#e38191,#cc607d,#ad466c,#8b3058,#672044
  
  #crispra
  #darkseagreen-5
  tdata_moi %>% 
    filter(system == sys_name) %>%
    ggplot(aes(x = volume_ul, y = moi, color = target_gene, group = target_gene)) +
    annotate(
      "rect",
      xmin = 1e-5,
      xmax = 10,
      ymin = 0.2,
      ymax = 0.4,
      fill = "gray85",
      alpha = 0.5
    ) +
    geom_hline(
      yintercept = c(0.2, 0.4),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.5
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_x_log10(
      breaks = axis_breaks,
      labels = axis_labels
    ) +
    scale_color_manual(values = target_colors) +
    labs(
      x = "Lentiviral input (µL)",
      y = "Estimated MOI",
      color = "Guide target"
    ) +
    theme_nathan()
}

plot_a <- plot_moi("CRISPRa")
plot_i <- plot_moi("CRISPRi")

print(plot_a)
print(plot_i)

ggsave(
  "CRISPRa_First_MOI_Experiment.png",
  plot = plot_a,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "CRISPRa_First_MOI_Experiment.pdf",
  plot = plot_a,
  width = 4,
  height = 2.5
)

ggsave(
  "CRISPRi_First_MOI_Experiment.png",
  plot = plot_i,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "CRISPRi_First_MOI_Experiment.pdf",
  plot = plot_i,
  width = 4,
  height = 2.5
)












#analyzing second experiment data - eek
library(tidyverse)
library(readr)

moi_raw <- read_csv(
  "viral titration day7 all well plates - Sheet1 copy.csv",
  col_names = FALSE,
  col_types = cols(.default = col_character())
)

# Keep only data rows
moi_df <- moi_raw %>%
  filter(str_detect(X1, "CRISPR")) %>%
  select(
    condition = X1,
    target_gene = X2,
    vol_10 = X3,
    vol_1 = X4,
    vol_n1 = X5,
    vol_n2 = X6,
    vol_n3 = X7,
    vol_n4 = X8,
    vol_n5 = X9,
    no_puro = X10,
    no_virus = X11
  ) %>%
  mutate(
    across(
      vol_10:no_virus,
      ~ parse_number(.x, locale = locale(grouping_mark = ","))
    ),
    target_gene = recode(
      target_gene,
      "G" = "GFI1B",
      "C" = "CD46",
      "NTC" = "NTC"
    )
  ) %>%
  separate(condition, into = c("system", "replicate"), sep = " - ")

# Impute missing controls from plate-matched values
# This uses the median no_puro/no_virus within each system + replicate
moi_df <- moi_df %>%
  group_by(system, replicate) %>%
  mutate(
    no_puro_used = if_else(
      is.na(no_puro),
      median(no_puro, na.rm = TRUE),
      no_puro
    ),
    no_virus_used = if_else(
      is.na(no_virus),
      median(no_virus, na.rm = TRUE),
      no_virus
    )
  ) %>%
  ungroup()

vol_map <- c(
  "vol_10" = 10,
  "vol_1" = 1,
  "vol_n1" = 0.1,
  "vol_n2" = 0.01,
  "vol_n3" = 0.001,
  "vol_n4" = 0.0001,
  "vol_n5" = 0.00001
)

moi_long <- moi_df %>%
  pivot_longer(
    cols = starts_with("vol_"),
    names_to = "volume_label",
    values_to = "signal"
  ) %>%
  mutate(
    volume_ul = vol_map[volume_label],
    pct_survival = 100 * (signal - no_virus_used) / (no_puro_used - no_virus_used),
    pct_survival = pmax(0, pmin(100, pct_survival)),
    moi = -log(1 - pmin(pct_survival / 100, 0.999))
  )

# Summarize across replicates
moi_summary <- moi_long %>%
  group_by(system, target_gene, volume_ul) %>%
  summarise(
    mean_survival = mean(pct_survival, na.rm = TRUE),
    se_survival = sd(pct_survival, na.rm = TRUE) / sqrt(sum(!is.na(pct_survival))),
    mean_moi = mean(moi, na.rm = TRUE),
    se_moi = sd(moi, na.rm = TRUE) / sqrt(sum(!is.na(moi))),
    n = sum(!is.na(moi)),
    .groups = "drop"
  )

# Find points in or near target range
target_points <- moi_summary %>%
  mutate(distance_from_0.3 = abs(mean_moi - 0.3)) %>%
  arrange(system, target_gene, distance_from_0.3) %>%
  group_by(system, target_gene) %>%
  slice(1) %>%
  ungroup()

print(target_points)

#plotting
axis_breaks <- c(10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001)
axis_labels <- c("10", "1", "10^-1", "10^-2", "10^-3", "10^-4", "10^-5")

plot_moi_summary <- function(sys_name) {
  
  target_colors <- if (sys_name == "CRISPRa") {
    c(
      "NTC" = "darkseagreen2",
      "GFI1B" = "darkseagreen4",
      "CD46" = "darkseagreen3"
    )
  } else {
    c(
      "NTC" = "#f4a3a8",
      "GFI1B" = "#8b3058",
      "CD46" = "#cc607d"
    )
  }
  
  moi_summary %>%
    filter(system == sys_name) %>%
    ggplot(aes(x = volume_ul, y = mean_moi, color = target_gene, group = target_gene)) +
    annotate(
      "rect",
      xmin = 1e-5,
      xmax = 10,
      ymin = 0.2,
      ymax = 0.4,
      fill = "gray85",
      alpha = 0.5
    ) +
    geom_hline(
      yintercept = c(0.2, 0.4),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.5
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(
        ymin = mean_moi - se_moi,
        ymax = mean_moi + se_moi
      ),
      width = 0.05,
      linewidth = 0.4
    ) +
    scale_x_log10(
      breaks = axis_breaks,
      labels = axis_labels
    ) +
    scale_color_manual(values = target_colors) +
    labs(
      x = "Lentiviral input (µL)",
      y = "Estimated MOI",
      color = "Guide target"
    ) +
    theme_nathan()
}

plot_a <- plot_moi_summary("CRISPRa")
plot_i <- plot_moi_summary("CRISPRi")

print(plot_a)
print(plot_i)

ggsave(
  "CRISPRa_2nd_MOI_Experiment.png",
  plot = plot_a,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "CRISPRa_2nd_MOI_Experiment.pdf",
  plot = plot_a,
  width = 8,
  height = 5
)

ggsave(
  "CRISPRi_2nd_MOI_Experiment.png",
  plot = plot_i,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  "CRISPRi_2nd_MOI_Experiment.pdf",
  plot = plot_i,
  width = 8,
  height = 5
)









# ============================================================================
# Compare MOI plots: linear y-axis vs log y-axis
# ============================================================================

library(tidyverse)
library(patchwork)

# X-axis setup
axis_breaks <- c(10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001)

axis_labels <- c(
  "10",
  "1",
  expression(10^-1),
  expression(10^-2),
  expression(10^-3),
  expression(10^-4),
  expression(10^-5)
)

# Y-axis setup for log-scale MOI
y_axis_breaks <- c(0.001, 0.01, 0.1, 0.2, 0.4, 1, 3, 10)

y_axis_labels <- c(
  "0.001",
  "0.01",
  "0.1",
  "0.2",
  "0.4",
  "1",
  "3",
  "10"
)

# ============================================================================
# Linear y-axis version
# ============================================================================

plot_moi_summary_linear <- function(sys_name) {
  
  target_colors <- if (sys_name == "CRISPRa") {
    c(
      "NTC" = "darkseagreen2",
      "GFI1B" = "darkseagreen4",
      "CD46" = "darkseagreen3"
    )
  } else {
    c(
      "NTC" = "#f4a3a8",
      "GFI1B" = "#8b3058",
      "CD46" = "#cc607d"
    )
  }
  
  moi_summary %>%
    filter(system == sys_name) %>%
    ggplot(
      aes(
        x = volume_ul,
        y = mean_moi,
        color = target_gene,
        group = target_gene
      )
    ) +
    annotate(
      "rect",
      xmin = 1e-5,
      xmax = 10,
      ymin = 0.2,
      ymax = 0.4,
      fill = "gray85",
      alpha = 0.5
    ) +
    geom_hline(
      yintercept = c(0.2, 0.4),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.5
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(
        ymin = pmax(mean_moi - se_moi, 0),
        ymax = mean_moi + se_moi
      ),
      width = 0.05,
      linewidth = 0.4
    ) +
    scale_x_log10(
      breaks = axis_breaks,
      labels = axis_labels
    ) +
    scale_color_manual(values = target_colors) +
    labs(
      x = "Lentiviral input (µL)",
      y = "Estimated MOI",
      color = "Guide target"
    ) +
    theme_nathan()
}

# ============================================================================
# Log y-axis version
# Important: log axes cannot display 0, so a small plotting floor is used.
# This only affects visualization, not the calculated MOI values.
# ============================================================================

plot_moi_summary_log <- function(sys_name) {
  
  target_colors <- if (sys_name == "CRISPRa") {
    c(
      "NTC" = "darkseagreen2",
      "GFI1B" = "darkseagreen4",
      "CD46" = "darkseagreen3"
    )
  } else {
    c(
      "NTC" = "#f4a3a8",
      "GFI1B" = "#8b3058",
      "CD46" = "#cc607d"
    )
  }
  
  moi_summary %>%
    filter(system == sys_name) %>%
    mutate(
      mean_moi_plot = pmax(mean_moi, 0.001),
      ymin_moi_plot = pmax(mean_moi - se_moi, 0.001),
      ymax_moi_plot = pmax(mean_moi + se_moi, 0.001)
    ) %>%
    ggplot(
      aes(
        x = volume_ul,
        y = mean_moi_plot,
        color = target_gene,
        group = target_gene
      )
    ) +
    annotate(
      "rect",
      xmin = 1e-5,
      xmax = 10,
      ymin = 0.2,
      ymax = 0.4,
      fill = "gray85",
      alpha = 0.5
    ) +
    geom_hline(
      yintercept = c(0.2, 0.4),
      linetype = "dashed",
      color = "gray40",
      linewidth = 0.5
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(
        ymin = ymin_moi_plot,
        ymax = ymax_moi_plot
      ),
      width = 0.05,
      linewidth = 0.4
    ) +
    scale_x_log10(
      breaks = axis_breaks,
      labels = axis_labels
    ) +
    scale_y_log10(
      breaks = y_axis_breaks,
      labels = y_axis_labels
    ) +
    scale_color_manual(values = target_colors) +
    labs(
      x = "Lentiviral input (µL)",
      y = "Estimated MOI (log scale)",
      color = "Guide target"
    ) +
    theme_nathan()
}

# ============================================================================
# Generate comparison plots
# ============================================================================

plot_a_linear <- plot_moi_summary_linear("CRISPRa")
plot_a_log <- plot_moi_summary_log("CRISPRa")

plot_i_linear <- plot_moi_summary_linear("CRISPRi")
plot_i_log <- plot_moi_summary_log("CRISPRi")

# Side-by-side comparison
moi_crispra_compare <- plot_a_linear + plot_a_log
moi_crispri_compare <- plot_i_linear + plot_i_log

print(moi_crispra_compare)
print(moi_crispri_compare)

# ============================================================================
# Save comparison plots
# ============================================================================

ggsave(
  "CRISPRa_MOI_og_yaxis.png",
  plot = plot_a_log,
  width = 12,
  height = 5,
  dpi = 300
)

ggsave(
  "CRISPRa_MOI_log_yaxis.pdf",
  plot = plot_a_log,
  width = 12,
  height = 5
)

ggsave(
  "CRISPRi_MOI_linear_yaxis.png",
  plot = plot_i_log,
  width = 12,
  height = 5,
  dpi = 300
)

ggsave(
  "CRISPRi_MOI_linear_yaxis.pdf",
  plot = plot_i_log,
  width = 12,
  height = 5
)
