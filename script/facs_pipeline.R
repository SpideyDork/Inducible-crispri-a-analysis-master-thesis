# ============================================================================
# FACS Pipeline — FCS gating, QC, and CD46/tdTomato plotting
# Combines flow_analysis.R + Flow_cytometry_plotting.R
# ============================================================================

library(flowCore)
library(ggcyto)
library(ggridges)
library(tidyverse)
library(gridExtra)

source("nathan_functs.R")

output_dir <- "final_images"

# ── Colors & helpers ──────────────────────────────────────────

system_colors <- c("CRISPRi" = "#8b3058", "CRISPRa" = "darkseagreen4")

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

clean_log_labels <- function(x) {
  x <- format(x, scientific = FALSE, trim = TRUE)
  x <- sub("0+$", "", x)
  sub("\\.$", "", x)
}

cell_line_to_system <- function(x) {
  case_when(x == "iC" ~ "CRISPRi", x == "aC" ~ "CRISPRa", TRUE ~ NA_character_)
}

convert_dose_to_uM <- function(inducer, dose) {
  case_when(
    inducer == "Dox"  ~ dose / 1000,
    inducer == "TMP"  ~ dose,
    inducer == "CTRL" ~ 0,
    TRUE              ~ dose
  )
}

label_k <- function(x) paste0(x / 1000, "K")

# ============================================================================
# SECTION 1: LOAD FCS FILES
# ============================================================================

fcs_dir <- "~/Documents/Master Thesis/Raw Data and Code/FACS data"

fs <- read.flowSet(path = fcs_dir, pattern = "\\.fcs$", truncate_max_range = FALSE)

# ============================================================================
# SECTION 2: GATE DEFINITIONS & APPLICATION
# ============================================================================

cell_gate <- polygonGate(
  "FSC-A" = c(175000, 250000, 600000, 750000, 600000, 200000),
  "SSC-A" = c(50000,  20000,  40000,  120000, 200000, 150000),
  filterId = "cells"
)

singlet_gate <- polygonGate(
  "FSC-A" = c(150000, 125000, 750000, 775000),
  "FSC-H" = c(50000,  125000, 350000, 275000),
  filterId = "singlets"
)

tdtomato_gate <- rectangleGate("FL3-A" = c(4000, Inf), filterId = "tdTomato_pos")

# Polygon coords for ggplot overlays
cell_gate_df    <- data.frame(x = c(125000,200000,600000,750000,600000,175000),
                               y = c(50000, 20000, 40000, 120000,250000,200000))
singlet_gate_df <- data.frame(x = c(150000,125000,750000,775000),
                               y = c(50000, 125000,350000,275000))

fs_cells      <- Subset(fs, cell_gate)
fs_singlets   <- Subset(fs_cells, singlet_gate)
fs_tomato_pos <- Subset(fs_singlets, tdtomato_gate)

extract_fs <- function(flowset) {
  lapply(sampleNames(flowset), function(nm) {
    df <- as.data.frame(exprs(flowset[[nm]]))
    if (nrow(df) == 0) return(NULL)
    df$sample <- sub("_Data Source.*", "", nm)
    df
  }) |> bind_rows()
}

df_raw        <- extract_fs(fs)
df_cells      <- extract_fs(fs_cells)
df_singlets   <- extract_fs(fs_singlets)
df_tomato_pos <- extract_fs(fs_tomato_pos)

# ============================================================================
# SECTION 3: SAMPLE ORDERING
# ============================================================================

ic_ctrl <- c("iC CTRL R1", "iC CTRL 2", "iC CTRL - 3")
ac_ctrl <- c("aC Ctrl - 1", "aC Ctrl - 2", "aC Ctrl - 3")

levels_iC_TMP <- c(ic_ctrl,
  paste("iC T0.01", c("- 1","- 2","- 3")), paste("iC T0.1", c("- 1","- 2","- 3")),
  paste("iC T1",   c("- 1","- 2","- 3")), paste("iC T5",   c("- 1","- 2","- 3")),
  paste("iC T10",  c("- 1","- 2","- 3")))

levels_iC_Dox <- c(ic_ctrl,
  paste("iC D25",   c("- 1","- 2","- 3")), paste("iC D50",   c("- 1","- 2","- 3")),
  paste("iC D200",  c("- 1","- 2","- 3")), paste("iC D600",  c("- 1","- 2","- 3")),
  paste("iC D2000", c("- 1","- 2","- 3")))

levels_aC_Dox <- c(ac_ctrl,
  paste("aC D25",   c("- 1","- 2","- 3")), paste("aC D50",   c("- 1","- 2","- 3")),
  paste("aC D200",  c("- 1","- 2","- 3")), paste("aC D600",  c("- 1","- 2","- 3")),
  paste("aC D2000", c("- 1","- 2","- 3")))

experimental_samples <- c(levels_iC_TMP, levels_iC_Dox, levels_aC_Dox)

order_samples <- function(df, levels) {
  df |> filter(sample %in% levels) |> mutate(sample = factor(sample, levels = levels))
}

# ============================================================================
# SECTION 4: PER-SAMPLE STATS & NORMALISATION
# ============================================================================

# % tdTomato+ per sample
n_singlets <- tibble(sample = sampleNames(fs_singlets)) |>
  mutate(
    sample     = sub("_Data Source.*", "", sample),
    n_singlets = map_int(sampleNames(fs_singlets), ~ nrow(exprs(fs_singlets[[.x]])))
  )

n_tomato <- tibble(sample = sampleNames(fs_tomato_pos)) |>
  mutate(
    sample   = sub("_Data Source.*", "", sample),
    n_tomato = map_int(sampleNames(fs_tomato_pos), ~ nrow(exprs(fs_tomato_pos[[.x]])))
  )

pct_cas9 <- left_join(n_singlets, n_tomato, by = "sample") |>
  mutate(
    n_tomato  = replace_na(n_tomato, 0),
    pct_cas9  = (n_tomato / n_singlets) * 100,
    cell_line = case_when(str_detect(sample, "^iC") ~ "iC", str_detect(sample, "^aC") ~ "aC"),
    inducer   = case_when(
      str_detect(sample, "CTRL|Ctrl") ~ "CTRL",
      str_detect(sample, " T")        ~ "TMP",
      str_detect(sample, " D")        ~ "Dox"
    ),
    dose = case_when(
      str_detect(sample, "CTRL|Ctrl") ~ 0,
      str_detect(sample, " T")        ~ as.numeric(str_extract(sample, "(?<=T)[0-9.]+")),
      str_detect(sample, " D")        ~ as.numeric(str_extract(sample, "(?<=D)[0-9]+"))
    )
  ) |>
  filter(!is.na(cell_line))

pct_cas9_summary <- pct_cas9 |>
  group_by(cell_line, inducer, dose) |>
  summarise(mean_pct = mean(pct_cas9), se_pct = se(pct_cas9), .groups = "drop")

# APC mean per sample (in tdTomato+ cells)
per_sample_stats <- df_tomato_pos |>
  filter(sample %in% experimental_samples) |>
  group_by(sample) |>
  summarise(
    sample_mean_APC = mean(`FL4-A`, na.rm = TRUE),
    median_APC      = median(`FL4-A`, na.rm = TRUE),
    se_APC          = se(`FL4-A`),
    n_events        = n(),
    .groups = "drop"
  ) |>
  mutate(
    cell_line = case_when(str_detect(sample, "^iC") ~ "iC", str_detect(sample, "^aC") ~ "aC"),
    inducer   = case_when(
      str_detect(sample, "CTRL|Ctrl") ~ "CTRL",
      str_detect(sample, " T")        ~ "TMP",
      str_detect(sample, " D")        ~ "Dox"
    ),
    dose = case_when(
      str_detect(sample, "CTRL|Ctrl") ~ 0,
      str_detect(sample, " T")        ~ as.numeric(str_extract(sample, "(?<=T)[0-9.]+")),
      str_detect(sample, " D")        ~ as.numeric(str_extract(sample, "(?<=D)[0-9]+"))
    ),
    replicate = str_extract(sample, "[0-9]$")
  )

ctrl_mean_iC <- per_sample_stats |> filter(cell_line == "iC", dose == 0) |> pull(sample_mean_APC) |> mean(na.rm = TRUE)
ctrl_mean_aC <- per_sample_stats |> filter(cell_line == "aC", dose == 0) |> pull(sample_mean_APC) |> mean(na.rm = TRUE)

per_sample_stats <- per_sample_stats |>
  mutate(
    system    = cell_line_to_system(cell_line),
    dose_uM   = convert_dose_to_uM(inducer, dose),
    dose_plot = case_when(
      inducer == "Dox"  & dose == 0 ~ 0.0125,
      inducer == "TMP"  & dose == 0 ~ 0.003,
      TRUE ~ dose_uM
    ),
    log2fc_APC = case_when(
      cell_line == "iC" ~ log2(sample_mean_APC / ctrl_mean_iC),
      cell_line == "aC" ~ log2(sample_mean_APC / ctrl_mean_aC)
    )
  )

per_condition_stats <- per_sample_stats |>
  group_by(cell_line, system, inducer, dose, dose_uM, dose_plot) |>
  summarise(
    mean_APC    = mean(sample_mean_APC, na.rm = TRUE),
    se_APC      = se(sample_mean_APC),
    mean_log2fc = mean(log2fc_APC, na.rm = TRUE),
    se_log2fc   = se(log2fc_APC),
    .groups = "drop"
  ) |>
  arrange(cell_line, inducer, dose)

# ============================================================================
# SECTION 5: QC GATE PLOTS (call manually as needed)
# ============================================================================

plot_cell_gate <- function(df) {
  ggplot(df, aes(`FSC-A`, `SSC-A`)) +
    geom_hex(bins = 128) + scale_fill_viridis_c(option = "G") +
    geom_polygon(data = cell_gate_df, aes(x, y), fill = NA, color = "red",
                 linewidth = 0.4, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(0,1e6), labels = label_k) +
    scale_y_continuous(limits = c(0,1e6), labels = label_k) +
    facet_wrap(~sample, ncol = 3) +
    labs(x = "FSC-A", y = "SSC-A") + theme_nathan()
}

plot_singlet_gate <- function(df) {
  ggplot(df, aes(`FSC-A`, `FSC-H`)) +
    geom_hex(bins = 128) + scale_fill_viridis_c(option = "G") +
    geom_polygon(data = singlet_gate_df, aes(x, y), fill = NA, color = "red",
                 linewidth = 0.4, inherit.aes = FALSE) +
    scale_x_continuous(limits = c(0,1e6), labels = label_k) +
    scale_y_continuous(limits = c(0,1e6), labels = label_k) +
    facet_wrap(~sample, ncol = 3) +
    labs(x = "FSC-A", y = "FSC-H") + theme_nathan()
}

plot_tdtomato_gate <- function(df) {
  ggplot(df, aes(`FL3-A`, `SSC-A`)) +
    geom_hex(bins = 128) + scale_fill_viridis_c(option = "G") +
    geom_vline(xintercept = 4000, color = "red", linewidth = 0.4, linetype = "dashed") +
    scale_x_log10(limits = c(10,1e8), labels = label_k) +
    scale_y_continuous(limits = c(0,1e6), labels = label_k) +
    facet_wrap(~sample, ncol = 3) +
    labs(x = "tdTomato (FL3-A)", y = "SSC-A") + theme_nathan()
}

plot_apc <- function(df) {
  ggplot(df, aes(`FL4-A`, `SSC-A`)) +
    geom_hex(bins = 128) + scale_fill_viridis_c(option = "G") +
    scale_x_log10(limits = c(10, 20000), labels = label_k) +
    scale_y_continuous(limits = c(0,1e6), labels = label_k) +
    facet_wrap(~sample, ncol = 3) +
    labs(x = "CD46 APC (FL4-A)", y = "SSC-A") + theme_nathan()
}

# Example QC calls:
# plot_cell_gate(order_samples(df_raw, levels_iC_Dox))
# plot_singlet_gate(order_samples(df_cells, levels_iC_Dox))
# plot_tdtomato_gate(order_samples(df_singlets, levels_iC_Dox))
# plot_apc(order_samples(df_tomato_pos, levels_iC_Dox))

# ============================================================================
# SECTION 6: GATING PIPELINE FIGURE (single control sample)
# ============================================================================

ctrl_raw      <- as.data.frame(exprs(fs[["iC CTRL R1_Data Source - 1.fcs"]]))
ctrl_cells    <- as.data.frame(exprs(fs_cells[["iC CTRL R1_Data Source - 1.fcs"]]))
ctrl_singlets <- as.data.frame(exprs(fs_singlets[["iC CTRL R1_Data Source - 1.fcs"]]))
ctrl_tomato   <- as.data.frame(exprs(fs_tomato_pos[["iC CTRL R1_Data Source - 1.fcs"]]))

p_gate1 <- ggplot(ctrl_raw, aes(`FSC-A`, `SSC-A`)) +
  geom_hex(bins = 128) + scale_fill_viridis_c(option = "G", guide = "none") +
  geom_polygon(data = cell_gate_df, aes(x, y), fill = NA, color = "red",
               linewidth = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0,1e6), labels = label_k) +
  scale_y_continuous(limits = c(0,1e6), labels = label_k) +
  labs(title = "Gate 1: Live cells", x = "FSC-A", y = "SSC-A") +
  theme_nathan()

p_gate2 <- ggplot(ctrl_cells, aes(`FSC-A`, `FSC-H`)) +
  geom_hex(bins = 128) + scale_fill_viridis_c(option = "G", guide = "none") +
  geom_polygon(data = singlet_gate_df, aes(x, y), fill = NA, color = "red",
               linewidth = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(0,1e6), labels = label_k) +
  scale_y_continuous(limits = c(0,1e6), labels = label_k) +
  labs(title = "Gate 2: Singlets", x = "FSC-A", y = "FSC-H") +
  theme_nathan()

p_gate3 <- ggplot(ctrl_singlets, aes(`FL3-A`, `SSC-A`)) +
  geom_hex(bins = 128) + scale_fill_viridis_c(option = "G", guide = "none") +
  geom_vline(xintercept = 4000, color = "red", linewidth = 0.5, linetype = "dashed") +
  scale_x_log10(limits = c(10,1e8), labels = label_k) +
  scale_y_continuous(limits = c(0,1e6), labels = label_k) +
  labs(title = "Gate 3: tdTomato-positive", x = "tdTomato (FL3-A)", y = "SSC-A") +
  theme_nathan()

p_gate4 <- ggplot(ctrl_tomato, aes(`FL4-A`, `SSC-A`)) +
  geom_hex(bins = 128) + scale_fill_viridis_c(option = "G", guide = "none") +
  scale_x_log10(limits = c(10,1e8), labels = label_k) +
  scale_y_continuous(limits = c(0,1e6), labels = label_k) +
  labs(title = "Gate 4: CD46 APC measurement", x = "APC CD46 (FL4-A)", y = "SSC-A") +
  theme_nathan()

fig_gating_pipeline <- arrangeGrob(p_gate1, p_gate2, p_gate3, p_gate4, nrow = 2)

# ============================================================================
# SECTION 6.1: Singlet Datafraes to compare TdTomato
# ============================================================================

# Compute % each side of gate for each control
pct_label <- function(df) {
  total <- nrow(df)
  pct_pos <- round(sum(df$`FL3-A` > 4000) / total * 100, 1)
  pct_neg <- round(100 - pct_pos, 1)
  list(pos = pct_pos, neg = pct_neg)
}

# Create control singlet dataframes for tdTomato comparison

ctrl_iC_singlets <- as.data.frame(
  exprs(fs_singlets[["iC CTRL R1_Data Source - 1.fcs"]])
)

ctrl_aC_singlets <- as.data.frame(
  exprs(fs_singlets[["aC Ctrl - 1_Data Source - 1.fcs"]])
)

iC_pct <- pct_label(ctrl_iC_singlets)
aC_pct <- pct_label(ctrl_aC_singlets)

p_iC_tomato <- ggplot(ctrl_iC_singlets, aes(x = `FL3-A`, y = `SSC-A`)) +
  geom_hex(bins = 128) +
  scale_fill_viridis_c(option = "G", guide = "none") +
  geom_vline(xintercept = 4000, color = "red", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = 500,        y = 950000, label = paste0(iC_pct$neg, "%"), color = "red") +
  annotate("text", x = 10000000,   y = 950000, label = paste0(iC_pct$pos, "%"), color = "red") +
  scale_x_log10(limits = c(10, 100000000), labels = label_k) +
  scale_y_continuous(limits = c(0, 1000000), labels = label_k) +
  labs(
    #title = "CRISPRi Control",
    x = "tdTomato (FL3-A)",
    y = "SSC-A"
  ) +
  theme_nathan_grid() +
  theme (plot.title = element_text(size = 10, face = "bold", margin = margin(b = 8)))

p_aC_tomato <- ggplot(ctrl_aC_singlets, aes(x = `FL3-A`, y = `SSC-A`)) +
  geom_hex(bins = 128) +
  scale_fill_viridis_c(option = "G", guide = "none") +
  geom_vline(xintercept = 4000, color = "red", linewidth = 0.5, linetype = "dashed") +
  annotate("text", x = 500,        y = 950000, label = paste0(aC_pct$neg, "%"), color = "red") +
  annotate("text", x = 10000000,   y = 950000, label = paste0(aC_pct$pos, "%"), color = "red") +
  scale_x_log10(limits = c(10, 100000000), labels = label_k) +
  scale_y_continuous(limits = c(0, 1000000), labels = label_k) +
  labs(
    #title = "CRISPRa Control",
    x = "tdTomato (FL3-A)",
    y = "SSC-A"
  ) +
  theme_nathan_grid() +
  theme (plot.title = element_text(size = 10, face = "bold", margin = margin(b = 8)))

grid.arrange(
  p_iC_tomato,
  p_aC_tomato,
  nrow = 1
  #top = "tdTomato-Positive Fraction in CRISPRa Control Cell Lines"
)

fig_tomato_pos_percent <- arrangeGrob(
  p_iC_tomato,
  p_aC_tomato,
  nrow = 1
)

# ============================================================================
# SECTION 7: RIDGE PLOTS
# ============================================================================

rep1_names <- c("iC CTRL R1", "aC Ctrl - 1", str_subset(unique(df_tomato_pos$sample), "- 1$"))

df_ridges <- df_tomato_pos |>
  filter(sample %in% rep1_names) |>
  mutate(
    cell_line = case_when(str_detect(sample, "^iC") ~ "iC", str_detect(sample, "^aC") ~ "aC"),
    inducer   = case_when(
      str_detect(sample, "CTRL|Ctrl") ~ "CTRL",
      str_detect(sample, " T")        ~ "TMP",
      str_detect(sample, " D")        ~ "Dox"
    ),
    dose      = case_when(
      str_detect(sample, "CTRL|Ctrl") ~ 0,
      str_detect(sample, " T")        ~ as.numeric(str_extract(sample, "(?<=T)[0-9.]+")),
      str_detect(sample, " D")        ~ as.numeric(str_extract(sample, "(?<=D)[0-9]+"))
    ),
    dose_uM   = convert_dose_to_uM(inducer, dose),
    condition = if_else(inducer == "CTRL", "Control", clean_log_labels(dose_uM))
  )

df_ridges_means <- df_ridges |>
  filter(is.finite(`FL4-A`), `FL4-A` > 0) |>
  group_by(cell_line, inducer, condition, dose, dose_uM) |>
  summarise(mean_FL4 = mean(`FL4-A`, na.rm = TRUE), .groups = "drop")

make_ridge_plot <- function(cl, ind) {
  y_label <- case_when(
    cl == "iC" & ind == "Dox" ~ "CRISPRi Dox conditions (\u00b5M)",
    cl == "iC" & ind == "TMP" ~ "CRISPRi TMP conditions (\u00b5M)",
    cl == "aC" & ind == "Dox" ~ "CRISPRa Dox conditions (\u00b5M)"
  )
  samp_df <- df_ridges |>
    filter(cell_line == cl, inducer %in% c(ind,"CTRL"), is.finite(`FL4-A`), `FL4-A` > 0)
  mean_df <- df_ridges_means |>
    filter(cell_line == cl, inducer %in% c(ind,"CTRL"), is.finite(mean_FL4), mean_FL4 > 0)
  cond_levels <- samp_df |> arrange(dose) |> pull(condition) |> unique()
  samp_df <- mutate(samp_df, condition = factor(condition, levels = cond_levels))
  mean_df <- mutate(mean_df, condition = factor(condition, levels = cond_levels))

  ggplot(samp_df, aes(`FL4-A`, condition)) +
    geom_density_ridges(aes(fill = condition), jittered_points = TRUE,
                        position = position_points_jitter(width = 0.05, height = 0.15),
                        point_shape = 16, point_size = 0.5, point_alpha = 0.2,
                        alpha = 0.6, scale = 0.9, show.legend = FALSE) +
    geom_point(data = mean_df, aes(x = mean_FL4, y = condition, fill = condition),
               shape = 23, size = 3, color = "black",
               position = position_nudge(y = -0.15), show.legend = FALSE) +
    scale_x_log10(limits = c(NA, 50000), labels = scales::comma) +
    scale_y_discrete(limits = rev) +
    labs(x = "APC signal (FL4-A)", y = y_label) +
    theme_nathan() + theme(legend.position = "none")
}

p_ridge_iC_dox <- make_ridge_plot("iC", "Dox")
p_ridge_iC_tmp <- make_ridge_plot("iC", "TMP")
p_ridge_aC_dox <- make_ridge_plot("aC", "Dox")

fig_ridges <- arrangeGrob(p_ridge_iC_tmp, p_ridge_iC_dox, p_ridge_aC_dox, ncol = 3)

# ============================================================================
# SECTION 8: TDTOMATO RETENTION PLOTS (bar)
# ============================================================================

# Helper to build retention data for one inducer
make_retention_data <- function(cl_vec, ind, dose_divisor = 1, ctrl_dose_plot = 0.0125) {
  dat <- bind_rows(lapply(cl_vec, function(cl) {
    pct_cas9 |>
      filter(cell_line == cl, inducer %in% c("CTRL", ind)) |>
      mutate(system = cell_line_to_system(cl))
  })) |>
    mutate(dose_uM = dose / dose_divisor)

  ctrl_ref <- dat |>
    filter(dose == 0) |>
    group_by(system) |>
    summarise(ctrl_mean_pct = mean(pct_cas9, na.rm = TRUE), .groups = "drop")

  dat |>
    left_join(ctrl_ref, by = "system") |>
    mutate(
      pct_of_control = (pct_cas9 / ctrl_mean_pct) * 100,
      system = factor(system, levels = c("CRISPRi","CRISPRa"))
    )
}

# Dox retention (Dox stored as nM → convert to µM)
dox_dose_labels <- c("0"="ctrl","25"="0.025","50"="0.05","200"="0.2","600"="0.6","2000"="2")

dox_retention_norm <- make_retention_data(c("iC","aC"), "Dox", dose_divisor = 1000) |>
  mutate(dose_label = factor(dox_dose_labels[as.character(dose)], levels = unname(dox_dose_labels)))

dox_retention_summary <- dox_retention_norm |>
  group_by(system, dose, dose_label) |>
  summarise(mean_pct_of_control = mean(pct_of_control, na.rm = TRUE),
            se_pct_of_control   = se(pct_of_control), .groups = "drop")

p_retention_dox <- ggplot() +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey60", linewidth = 0.7) +
  geom_col(data = dox_retention_summary,
           aes(dose_label, mean_pct_of_control, fill = system),
           position = position_dodge(width = 0.7), width = 0.65, alpha = 0.7) +
  geom_errorbar(data = dox_retention_summary,
                aes(dose_label, ymin = mean_pct_of_control - se_pct_of_control,
                    ymax = mean_pct_of_control + se_pct_of_control, group = system),
                position = position_dodge(width = 0.7), width = 0.2,
                linewidth = 0.6, color = "black") +
  geom_point(data = dox_retention_norm,
             aes(dose_label, pct_of_control, color = system),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7, seed = 42),
             size = 1.9, alpha = 0.75) +
  scale_fill_manual(values = system_colors) +
  scale_color_manual(values = system_colors) +
  scale_y_continuous(limits = c(0, 120), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Dox concentration (\u00b5M)", y = "tdTomato-positive cells (% of control)",
       fill = NULL, color = NULL) +
  theme_nathan() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# TMP retention (TMP already in µM)
tmp_dose_labels <- c("0"="ctrl","0.01"="0.01","0.1"="0.1","1"="1","5"="5","10"="10")

tmp_retention_norm <- make_retention_data("iC", "TMP", dose_divisor = 1) |>
  mutate(dose_label = factor(tmp_dose_labels[as.character(dose)], levels = unname(tmp_dose_labels)))

tmp_retention_summary <- tmp_retention_norm |>
  group_by(system, dose, dose_label) |>
  summarise(mean_pct_of_control = mean(pct_of_control, na.rm = TRUE),
            se_pct_of_control   = se(pct_of_control), .groups = "drop")

p_retention_tmp <- ggplot() +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey60", linewidth = 0.7) +
  geom_col(data = tmp_retention_summary,
           aes(dose_label, mean_pct_of_control, fill = system),
           width = 0.55, alpha = 0.7) +
  geom_errorbar(data = tmp_retention_summary,
                aes(dose_label, ymin = mean_pct_of_control - se_pct_of_control,
                    ymax = mean_pct_of_control + se_pct_of_control),
                width = 0.2, linewidth = 0.6, color = "black") +
  geom_point(data = tmp_retention_norm,
             aes(dose_label, pct_of_control, color = system),
             position = position_jitter(width = 0.1, height = 0, seed = 42),
             size = 1.9, alpha = 0.75) +
  scale_fill_manual(values = system_colors) +
  scale_color_manual(values = system_colors) +
  scale_y_continuous(limits = c(0, 120), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "TMP concentration (\u00b5M)", y = "tdTomato-positive cells (% of control)",
       fill = NULL, color = NULL) +
  theme_nathan() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================================
# SECTION 9: CD46 FOLD CHANGE PLOTS (FACS-based)
# ============================================================================

compare_df <- per_sample_stats |>
  left_join(pct_cas9 |> select(sample, pct_cas9), by = "sample")

make_cd46_fold_change_plot <- function(cl_vec, ind, ctrl_dose_plot, dose_divisor = 1) {
  dat <- bind_rows(lapply(cl_vec, function(cl) {
    compare_df |>
      filter(cell_line == cl, inducer %in% c("CTRL", ind)) |>
      mutate(system = cell_line_to_system(cl))
  })) |>
    mutate(
      dose_uM   = dose / dose_divisor,
      dose_plot = if_else(dose == 0, ctrl_dose_plot, dose_uM),
      system    = factor(system, levels = c("CRISPRi","CRISPRa"))
    )

  summary_dat <- dat |>
    group_by(system, dose, dose_plot) |>
    summarise(mean_log2fc = mean(log2fc_APC, na.rm = TRUE),
              se_log2fc   = se(log2fc_APC), .groups = "drop")

  x_breaks <- sort(unique(summary_dat$dose_plot))
  x_labels <- clean_log_labels(x_breaks)
  x_labels[x_breaks == ctrl_dose_plot] <- "ctrl"
  x_lab <- if (ind == "Dox") "Dox concentration (\u00b5M, log scale)" else "TMP concentration (\u00b5M, log scale)"

  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey65") +
    geom_line(data = summary_dat, aes(dose_plot, mean_log2fc, color = system, group = system),
              linewidth = 0.8) +
    geom_errorbar(data = summary_dat,
                  aes(dose_plot, ymin = mean_log2fc - se_log2fc, ymax = mean_log2fc + se_log2fc),
                  width = if (ind == "Dox") 0.015 else 0.04, linewidth = 0.6, color = "black") +
    geom_point(data = summary_dat, aes(dose_plot, mean_log2fc, color = system), size = 2.8) +
    scale_x_log10(breaks = x_breaks, labels = x_labels) +
    scale_color_manual(values = system_colors) +
    labs(x = x_lab, y = bquote(CD46 ~ log[2] ~ fold ~ change), color = NULL) +
    theme_nathan() +
    theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1), legend.position = "none")
}

p_dox_cd46 <- make_cd46_fold_change_plot(c("iC","aC"), "Dox", ctrl_dose_plot = 0.0125, dose_divisor = 1000)
p_tmp_cd46 <- make_cd46_fold_change_plot("iC", "TMP", ctrl_dose_plot = 0.003, dose_divisor = 1)

# ============================================================================
# SECTION 10: SAVE
# ============================================================================

plots_individual <- list(
  p_retention_dox = p_retention_dox,
  p_retention_tmp = p_retention_tmp,
  p_dox_cd46      = p_dox_cd46,
  p_tmp_cd46      = p_tmp_cd46,
  p_ridge_iC_dox  = p_ridge_iC_dox,
  p_ridge_iC_tmp  = p_ridge_iC_tmp,
  p_ridge_aC_dox  = p_ridge_aC_dox
)

walk2(names(plots_individual), plots_individual, function(name, plot) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot, width = 4, height = 2.5, dpi = 300)
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot, width = 4, height = 2.5)
  cat("Saved:", name, "\n")
})

# Combined figures
ggsave(file.path(output_dir, "fig_gating_pipeline.png"), fig_gating_pipeline, width = 4, height = 2.5, dpi = 300)
ggsave(file.path(output_dir, "fig_gating_pipeline.pdf"), fig_gating_pipeline, width = 4, height = 2.5)
ggsave(file.path(output_dir, "tomato_pos_percent.png"), fig_tomato_pos_percent, width = 4, height = 2.5, dpi = 300)
ggsave(file.path(output_dir, "tomato_pos_percent.pdf"), fig_tomato_pos_percent, width = 4, height = 2.5)
ggsave(file.path(output_dir, "fig_ridges_all.png"),      fig_ridges,          width = 8, height = 4, dpi = 300)
ggsave(file.path(output_dir, "fig_ridges_all.pdf"),      fig_ridges,          width = 8, height = 4)
cat("Saved combined figures\n")

