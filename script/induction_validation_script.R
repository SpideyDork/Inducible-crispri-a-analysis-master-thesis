# ============================================================================
# Combined qPCR Analysis — CD46 & GFI1B Induction Validation
# Both-HK normalisation, shared y-axis
# ============================================================================

library(tidyverse)
source("nathan_functs.R")

output_dir <- "final_images"

# ── Colors ───────────────────────────────────────────────────

perturbation_colors <- c("CRISPRi" = "#8b3058", "CRISPRa" = "darkseagreen4")

# ── Helpers ───────────────────────────────────────────────────

calculate_cq <- function(amp_raw, threshold = 5000) {
  amp_raw %>%
    select(-1) %>%
    pivot_longer(cols = -Cycle, names_to = "Well", values_to = "Fluor") %>%
    mutate(
      Row = str_extract(Well, "^[A-H]"),
      Col = as.integer(str_extract(Well, "\\d+$"))
    ) %>%
    group_by(Well, Row, Col) %>%
    arrange(Cycle, .by_group = TRUE) %>%
    mutate(
      baseline = mean(Fluor[Cycle <= 10], na.rm = TRUE),
      dRn      = Fluor - baseline
    ) %>%
    summarise(
      Cq = {
        idx <- which(dRn >= threshold)[1]
        if (is.na(idx))       NA_real_
        else if (idx == 1)    Cycle[1]
        else {
          x0 <- Cycle[idx-1]; x1 <- Cycle[idx]
          y0 <- dRn[idx-1];   y1 <- dRn[idx]
          x0 + (threshold - y0) / (y1 - y0)
        }
      },
      .groups = "drop"
    )
}

compute_fold_change <- function(delta_cq) {
  ctrl_stats <- delta_cq %>%
    filter(Treatment == "Control") %>%
    group_by(Perturbation) %>%
    summarise(ctrl_mean = mean(deltaCq_both, na.rm = TRUE), .groups = "drop")

  delta_cq %>%
    filter(Treatment == "Treated") %>%
    left_join(ctrl_stats, by = "Perturbation") %>%
    mutate(
      log2FC      = -(deltaCq_both - ctrl_mean),
      fold_change = 2^log2FC
    )
}

# ── CD46 pipeline ────────────────────────────────────────────

amp_cd46 <- read_csv(
  "admin_2026-04-21 10-45-17_BR008471 -  Quantification Amplification Results_SYBR.csv",
  show_col_types = FALSE
)

cq_cd46 <- calculate_cq(amp_cd46) %>%
  filter(Row %in% LETTERS[1:6], Col %in% 1:12, !(Row == "F" & Col %in% 7:12)) %>%
  mutate(
    Sample = case_when(
      Row == "A" & Col %in% 1:6  ~ "Inhibition_Treated_R1",
      Row == "A" & Col %in% 7:12 ~ "Activation_Treated_R1",
      Row == "B" & Col %in% 1:6  ~ "Inhibition_Treated_R2",
      Row == "B" & Col %in% 7:12 ~ "Activation_Treated_R2",
      Row == "C" & Col %in% 1:6  ~ "Inhibition_Treated_R3",
      Row == "C" & Col %in% 7:12 ~ "Activation_Treated_R3",
      Row == "D" & Col %in% 1:6  ~ "Inhibition_Control_R1",
      Row == "D" & Col %in% 7:12 ~ "Activation_Control_R1",
      Row == "E" & Col %in% 1:6  ~ "Inhibition_Control_R2",
      Row == "E" & Col %in% 7:12 ~ "Activation_Control_R2",
      Row == "F" & Col %in% 1:6  ~ "Inhibition_Control_R3"
    ),
    Gene         = case_when(
      Col %in% c(1,2,7,8)   ~ "CD46",
      Col %in% c(3,4,9,10)  ~ "Actin",
      Col %in% c(5,6,11,12) ~ "18S"
    ),
    Perturbation = if_else(str_detect(Sample, "Inhibition"), "CRISPRi", "CRISPRa"),
    Treatment    = if_else(str_detect(Sample, "Control"), "Control", "Treated")
  ) %>%
  filter(Well != "B10")

# Drop tech replicates with Cq range >= 1, then average
bad_tech <- cq_cd46 %>%
  filter(!is.na(Cq)) %>%
  group_by(Sample, Gene) %>%
  summarise(range_Cq = max(Cq) - min(Cq), .groups = "drop") %>%
  filter(range_Cq >= 1)

mean_cq_cd46 <- cq_cd46 %>%
  filter(!is.na(Cq), !is.na(Sample), !is.na(Gene)) %>%
  anti_join(bad_tech, by = c("Sample", "Gene")) %>%
  group_by(Sample, Gene, Perturbation, Treatment) %>%
  summarise(mean_Cq = mean(Cq, na.rm = TRUE), .groups = "drop")

# Remove samples with failed HK genes
failed_hk <- mean_cq_cd46 %>%
  filter(Gene %in% c("Actin", "18S"), mean_Cq > 30) %>%
  pull(Sample)

fc_cd46 <- mean_cq_cd46 %>%
  filter(!Sample %in% failed_hk) %>%
  mutate(Bio_Rep = str_extract(Sample, "R(\\d+)$", group = 1)) %>%
  pivot_wider(names_from = Gene, values_from = mean_Cq) %>%
  mutate(deltaCq_both = CD46 - (Actin + `18S`) / 2) %>%
  compute_fold_change() %>%
  mutate(
    Target       = "CD46",
    Perturbation = factor(Perturbation, levels = c("CRISPRi", "CRISPRa"))
  ) %>%
  filter(!is.na(log2FC), is.finite(log2FC))

# ── GFI1B pipeline ───────────────────────────────────────────

amp_gfi1b <- read_csv(
  "admin_2026-04-07 12-10-51_BR008471 -  Quantification Amplification Results_SYBR.csv",
  show_col_types = FALSE
)

fc_gfi1b <- calculate_cq(amp_gfi1b) %>%
  filter(Row %in% LETTERS[1:6]) %>%
  mutate(
    Sample = case_when(
      Row == "A" & Col %in% 1:6  ~ "IGw1",
      Row == "A" & Col %in% 7:12 ~ "NTC2",
      Row == "B" & Col %in% 1:6  ~ "IGw2",
      Row == "B" & Col %in% 7:12 ~ "aGw1",
      Row == "C" & Col %in% 1:6  ~ "IGw3",
      Row == "C" & Col %in% 7:12 ~ "aGw2",
      Row == "D" & Col %in% 1:6  ~ "IGwo1",
      Row == "D" & Col %in% 7:12 ~ "aGwo1",
      Row == "E" & Col %in% 1:6  ~ "IGwo2",
      Row == "E" & Col %in% 7:12 ~ "aGwo2",
      Row == "F" & Col %in% 1:6  ~ "IGwo3",
      Row == "F" & Col %in% 7:12 ~ "aGwo3"
    ),
    Gene         = case_when(
      Col %in% c(1,2,7,8)   ~ "GFI1B",
      Col %in% c(3,4,9,10)  ~ "Actin",
      Col %in% c(5,6,11,12) ~ "18S"
    ),
    Perturbation = case_when(
      str_detect(Sample, "^IG") ~ "CRISPRi",
      str_detect(Sample, "^aG") ~ "CRISPRa"
    ),
    Treatment    = if_else(str_detect(Sample, "wo"), "Control", "Treated"),
    Bio_Rep      = str_extract(Sample, "\\d+$")
  ) %>%
  filter(!is.na(Perturbation), !is.na(Cq), !is.na(Gene)) %>%
  group_by(Sample, Gene, Perturbation, Treatment, Bio_Rep) %>%
  summarise(mean_Cq = mean(Cq, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Gene, values_from = mean_Cq) %>%
  mutate(deltaCq_both = GFI1B - (Actin + `18S`) / 2) %>%
  compute_fold_change() %>%
  mutate(
    Target       = "GFI1B",
    Perturbation = factor(Perturbation, levels = c("CRISPRi", "CRISPRa"))
  ) %>%
  filter(!is.na(log2FC), is.finite(log2FC))

# ── Shared y-axis ─────────────────────────────────────────────

all_fc   <- c(fc_cd46$log2FC, fc_gfi1b$log2FC)
y_margin <- diff(range(all_fc)) * 0.1
y_limits <- range(all_fc) + c(-y_margin, y_margin)

# ── Plot function ─────────────────────────────────────────────

make_plot <- function(fc_data, target, y_lim) {
  y_lab <- bquote(.(target) ~ log[2] ~ fold ~ change)

  summary_data <- fc_data %>%
    group_by(Perturbation) %>%
    summarise(
      mean_log2FC = mean(log2FC, na.rm = TRUE),
      se_log2FC   = sd(log2FC, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  ggplot() +
    geom_col(
      data = summary_data,
      aes(Perturbation, mean_log2FC, fill = Perturbation),
      width = 0.6, alpha = 0.7
    ) +
    geom_errorbar(
      data = summary_data,
      aes(Perturbation, ymin = mean_log2FC - se_log2FC, ymax = mean_log2FC + se_log2FC),
      width = 0.2, linewidth = 0.8, color = "black"
    ) +
    geom_point(
      data = fc_data,
      aes(Perturbation, log2FC, color = Perturbation, shape = Bio_Rep),
      position = position_jitter(width = 0.1, height = 0, seed = 42),
      size = 3, alpha = 0.85
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.6) +
    scale_fill_manual(values = perturbation_colors, guide = "none") +
    scale_color_manual(values = perturbation_colors, guide = "none") +
    scale_y_continuous(limits = y_lim) +
    labs(x = NULL, y = y_lab, shape = "Biological replicate") +
    theme_nathan() +
    theme(legend.position = "right")
}

p_cd46  <- make_plot(fc_cd46,  "CD46",  y_limits)
p_gfi1b <- make_plot(fc_gfi1b, "GFI1B", y_limits)

# ── Save ──────────────────────────────────────────────────────

plots <- list(
  CD46_induction_validation  = p_cd46,
  GFI1B_induction_validation = p_gfi1b
)

for (name in names(plots)) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plots[[name]], width = 4, height = 4, dpi = 300)
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plots[[name]], width = 4, height = 4)
  cat("Saved:", name, "\n")
}
