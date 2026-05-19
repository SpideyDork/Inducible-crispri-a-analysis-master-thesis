# ============================================================================
# Full qPCR pipeline — 14 plates → 6 thesis plots
# Outputs: 3 CD46 Hill fit plots + 3 GFI1B bar plots
# ============================================================================

library(tidyverse)
library(drc)

output_dir <- "final_images"

# ── Colors & theme ───────────────────────────────────────────

system_colors <- c("CRISPRi" = "#8b3058", "CRISPRa" = "darkseagreen4")

if (!exists("theme_nathan")) {
  theme_nathan <- function(base_size = 14) {
    theme_classic(base_size = base_size) +
      theme(
        axis.line        = element_line(linewidth = 0.6),
        axis.ticks       = element_line(linewidth = 0.6),
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold"),
        legend.title     = element_text(face = "plain"),
        legend.position  = "right"
      )
  }
}

# ============================================================================
# SECTION 1: PLATE FILES & DESIGN
# ============================================================================

default_threshold <- 5000

plate_files <- tribble(
  ~Plate,    ~file,                                                                                                                   ~Threshold,
  "plate1",  "admin_2026-05-02 12-43-26_BR008471_cd46_crispr1_plate1try2 -  Quantification Amplification Results_SYBR.csv",          default_threshold,
  "plate2",  "admin_2026-05-02 19-26-45_BR008471_indCRISPRi_cd46_plate2new -  Quantification Amplification Results_SYBR.csv",        default_threshold,
  "plate3",  "admin_2026-05-03 12-16-46_BR008471_indcrispri_cd46_plate3new -  Quantification Amplification Results_SYBR.csv",        default_threshold,
  "plate4",  "admin_2026-05-03 13-52-36_BR008471_indcrispri_cd46_plate4 -  Quantification Amplification Results_SYBR.csv",           default_threshold,
  "plate5",  "admin_2026-05-03 16-28-28_BR008471_indCRISPRi_cd46_plate5 -  Quantification Amplification Results_SYBR.csv",           default_threshold,
  "plate6",  "admin_2026-05-06 10-05-56_BR008471_indcrispri_gfi1b_plate6 -  Quantification Amplification Results_SYBR.csv",          default_threshold,
  "plate7",  "admin_2026-05-06 11-23-32_BR008471_indcrispri_gfi1b_plate7_t0.01_t10 -  Quantification Amplification Results_SYBR.csv",default_threshold,
  "plate8",  "admin_2026-05-04 09-45-58_BR008471_indCRISPRi_NTC_CD46_plate8 -  Quantification Amplification Results_SYBR.csv",       default_threshold,
  "plate9",  "admin_2026-05-04 12-12-18_BR008471_indCRISPRi_cd46_plate9_NTC -  Quantification Amplification Results_SYBR.csv",       default_threshold,
  "plate10", "admin_2026-05-04 15-09-00_BR008471_indCRISPRa_cd46_plate10_d25_d50 -  Quantification Amplification Results_SYBR.csv",  default_threshold,
  "plate11", "admin_2026-05-04 17-34-43_BR008471_indcrispra_cd46_plate11_d200_d600 -  Quantification Amplification Results_SYBR.csv",default_threshold,
  "plate12", "admin_2026-05-06 14-27-45_BR008471_indcrispra_cd46_gfi1b_plate12 -  Quantification Amplification Results_SYBR.csv",    default_threshold,
  "plate13", "admin_2026-05-06 12-37-09_BR008471_indcrispra_gfi1b_plate13 -  Quantification Amplification Results_SYBR.csv",         default_threshold,
  "plate14", "admin_2026-05-07 09-52-10_BR008471-indcrispri_a_cd46_gfi1b_plate14_ctrls -  Quantification Amplification Results_SYBR.csv", default_threshold
)

plate_3slot_design <- tribble(
  ~Plate,    ~Slot,    ~System,    ~Guide,   ~TargetGene, ~Condition,   ~DoseLabel,
  # CRISPRi CD46
  "plate1",  "slot1",  "CRISPRi",  "CD46",  "CD46",      "iC ctrl",   "iC ctrl",
  "plate1",  "slot2",  "CRISPRi",  "CD46",  "CD46",      "D25",       "D25",
  "plate1",  "slot3",  "CRISPRi",  "CD46",  "CD46",      "D50",       "D50",
  "plate2",  "slot1",  "CRISPRi",  "CD46",  "CD46",      "iC ctrl",   "iC ctrl",
  "plate2",  "slot2",  "CRISPRi",  "CD46",  "CD46",      "D200",      "D200",
  "plate2",  "slot3",  "CRISPRi",  "CD46",  "CD46",      "D600",      "D600",
  "plate3",  "slot1",  "CRISPRi",  "CD46",  "CD46",      "iC ctrl",   "iC ctrl",
  "plate3",  "slot2",  "CRISPRi",  "CD46",  "CD46",      "D2000",     "D2000",
  "plate3",  "slot3",  "CRISPRi",  "CD46",  "CD46",      "T0.01",     "T0.01",
  "plate4",  "slot1",  "CRISPRi",  "CD46",  "CD46",      "iC ctrl",   "iC ctrl",
  "plate4",  "slot2",  "CRISPRi",  "CD46",  "CD46",      "T0.1",      "T0.1",
  "plate4",  "slot3",  "CRISPRi",  "CD46",  "CD46",      "T1",        "T1",
  "plate5",  "slot1",  "CRISPRi",  "CD46",  "CD46",      "iC ctrl",   "iC ctrl",
  "plate5",  "slot2",  "CRISPRi",  "CD46",  "CD46",      "T5",        "T5",
  "plate5",  "slot3",  "CRISPRi",  "CD46",  "CD46",      "T10",       "T10",
  # CRISPRi GFI1B
  "plate6",  "slot1",  "CRISPRi",  "GFI1B", "GFI1B",     "iG ctrl",   "iG ctrl",
  "plate6",  "slot2",  "CRISPRi",  "GFI1B", "GFI1B",     "D50",       "D50",
  "plate6",  "slot3",  "CRISPRi",  "GFI1B", "GFI1B",     "D600",      "D600",
  "plate7",  "slot1",  "CRISPRi",  "GFI1B", "GFI1B",     "iG ctrl",   "iG ctrl",
  "plate7",  "slot2",  "CRISPRi",  "GFI1B", "GFI1B",     "T0.01",     "T0.01",
  "plate7",  "slot3",  "CRISPRi",  "GFI1B", "GFI1B",     "T10",       "T10",
  # CRISPRi NTC / toxicity
  "plate8",  "slot1",  "CRISPRi",  "CD46",  "CD46",      "iC ctrl",   "iC ctrl",
  "plate8",  "slot2",  "CRISPRi",  "NTC",   "CD46",      "iNTC ctrl", "iNTC ctrl",
  "plate8",  "slot3",  "CRISPRi",  "NTC",   "CD46",      "tox1",      "tox1",
  "plate9",  "slot1",  "CRISPRi",  "NTC",   "CD46",      "iNTC ctrl", "iNTC ctrl",
  "plate9",  "slot2",  "CRISPRi",  "NTC",   "CD46",      "tox2",      "tox2",
  "plate9",  "slot3",  "WT",       "WT",    "CD46",      "WT K562",   "WT K562",
  # CRISPRa CD46
  "plate10", "slot1",  "CRISPRa",  "CD46",  "CD46",      "aC ctrl",   "aC ctrl",
  "plate10", "slot2",  "CRISPRa",  "CD46",  "CD46",      "D25",       "D25",
  "plate10", "slot3",  "CRISPRa",  "CD46",  "CD46",      "D50",       "D50",
  "plate11", "slot1",  "CRISPRa",  "CD46",  "CD46",      "aC ctrl",   "aC ctrl",
  "plate11", "slot2",  "CRISPRa",  "CD46",  "CD46",      "D200",      "D200",
  "plate11", "slot3",  "CRISPRa",  "CD46",  "CD46",      "D600",      "D600",
  "plate12", "slot1",  "CRISPRa",  "CD46",  "CD46",      "aC ctrl",   "aC ctrl",
  "plate12", "slot2",  "CRISPRa",  "CD46",  "CD46",      "D2000",     "D2000",
  "plate12", "slot3",  "CRISPRa",  "GFI1B", "GFI1B",     "aG ctrl",   "aG ctrl",
  # CRISPRa GFI1B
  "plate13", "slot1",  "CRISPRa",  "GFI1B", "GFI1B",     "aG ctrl",   "aG ctrl",
  "plate13", "slot2",  "CRISPRa",  "GFI1B", "GFI1B",     "D50",       "D50",
  "plate13", "slot3",  "CRISPRa",  "GFI1B", "GFI1B",     "D600",      "D600",
  # Cross-plate controls
  "plate14", "slot1",  "CRISPRi",  "CD46",  "CD46",      "iC ctrl",   "iC ctrl",
  "plate14", "slot2",  "CRISPRa",  "CD46",  "CD46",      "aC ctrl",   "aC ctrl",
  "plate14", "slot3",  "CRISPRi",  "GFI1B", "GFI1B",     "iG ctrl",   "iG ctrl"
)

# ============================================================================
# SECTION 2: CQ EXTRACTION FUNCTIONS
# ============================================================================

read_qpcr_amp_long <- function(file, plate_id) {
  raw <- read_csv(file, show_col_types = FALSE)
  if (!"Cycle" %in% names(raw)) names(raw)[1] <- "Cycle"
  raw %>%
    select(Cycle, matches("^[A-H]([1-9]|1[0-2])$")) %>%
    pivot_longer(cols = -Cycle, names_to = "Well", values_to = "Fluor") %>%
    mutate(
      Plate = plate_id,
      Row   = str_extract(Well, "^[A-H]"),
      Col   = as.integer(str_extract(Well, "\\d+$"))
    )
}

make_3slot_map_one <- function(Plate, Slot, System, Guide, TargetGene, Condition, DoseLabel) {
  rows <- switch(Slot, slot1 = c("A","B","C"), slot2 = c("D","E","F"), slot3 = c("G","H"))

  regular_map <- tibble(Row = rows) %>%
    mutate(Bio_Rep = case_when(
      Row %in% c("A","D","G") ~ "R1",
      Row %in% c("B","E","H") ~ "R2",
      Row %in% c("C","F")     ~ "R3"
    )) %>%
    crossing(tibble(Col = 1:9)) %>%
    mutate(
      Gene     = case_when(Col %in% 1:3 ~ TargetGene, Col %in% 4:6 ~ "Actin", Col %in% 7:9 ~ "18S"),
      Tech_Rep = case_when(Col %in% c(1,4,7) ~ "T1", Col %in% c(2,5,8) ~ "T2", Col %in% c(3,6,9) ~ "T3")
    )

  special_map <- if (Slot == "slot3") {
    crossing(tibble(Row = c("A","B","C")), tibble(Col = 10:12)) %>%
      mutate(
        Bio_Rep  = "R3",
        Gene     = case_when(Row == "A" ~ TargetGene, Row == "B" ~ "Actin", Row == "C" ~ "18S"),
        Tech_Rep = case_when(Col == 10 ~ "T1", Col == 11 ~ "T2", Col == 12 ~ "T3")
      )
  } else tibble()

  bind_rows(regular_map, special_map) %>%
    mutate(
      Plate = Plate, Slot = Slot, System = System, Guide = Guide,
      TargetGene = TargetGene, Condition = Condition, DoseLabel = DoseLabel,
      Well   = paste0(Row, Col),
      Sample = paste(Plate, System, Guide, TargetGene, Condition, DoseLabel, Bio_Rep, sep = "__")
    ) %>%
    select(Plate, Well, Row, Col, Slot, System, Guide, TargetGene,
           Condition, DoseLabel, Bio_Rep, Tech_Rep, Gene, Sample)
}

add_analysis_metadata <- function(df) {
  df %>%
    mutate(
      CellLine = case_when(
        System == "CRISPRi" ~ "CRISPRi",
        System == "CRISPRa" ~ "CRISPRa",
        System == "WT"      ~ "WT"
      ),
      ControlGroup = case_when(
        Condition == "iNTC ctrl" ~ "iNTC", Condition == "iC ctrl" ~ "iC",
        Condition == "iG ctrl"   ~ "iG",   Condition == "aC ctrl" ~ "aC",
        Condition == "aG ctrl"   ~ "aG"
      ),
      IsControl = !is.na(ControlGroup),
      DoseType  = case_when(
        str_detect(DoseLabel, "^D")      ~ "Dox",
        str_detect(DoseLabel, "^T")      ~ "TMP",
        str_detect(DoseLabel, "ctrl")    ~ "Control",
        DoseLabel %in% c("tox1","tox2") ~ "Toxicity/NTC",
        DoseLabel == "WT K562"          ~ "WT",
        TRUE ~ "Other"
      ),
      DoseValue = case_when(
        str_detect(DoseLabel, "^D") ~ as.numeric(str_extract(DoseLabel, "(?<=^D)\\d+(\\.\\d+)?")),
        str_detect(DoseLabel, "^T") ~ as.numeric(str_extract(DoseLabel, "(?<=^T)\\d+(\\.\\d+)?")),
        str_detect(DoseLabel, "ctrl") ~ 0,
        TRUE ~ NA_real_
      ),
      NormalizeToCondition = case_when(
        Plate %in% c("plate8","plate9")        ~ "iNTC ctrl",
        System == "CRISPRi" & Guide == "CD46"  ~ "iC ctrl",
        System == "CRISPRi" & Guide == "GFI1B" ~ "iG ctrl",
        System == "CRISPRi" & Guide == "NTC"   ~ "iNTC ctrl",
        System == "CRISPRa" & Guide == "CD46"  ~ "aC ctrl",
        System == "CRISPRa" & Guide == "GFI1B" ~ "aG ctrl"
      ),
      NormalizeToGroup = case_when(
        NormalizeToCondition == "iNTC ctrl" ~ "iNTC",
        NormalizeToCondition == "iC ctrl"   ~ "iC",
        NormalizeToCondition == "iG ctrl"   ~ "iG",
        NormalizeToCondition == "aC ctrl"   ~ "aC",
        NormalizeToCondition == "aG ctrl"   ~ "aG"
      )
    )
}

calc_cq_one_plate <- function(Plate, file, Threshold) {
  read_qpcr_amp_long(file, Plate) %>%
    group_by(Plate, Well, Row, Col) %>%
    arrange(Cycle, .by_group = TRUE) %>%
    mutate(
      baseline = mean(Fluor[Cycle <= 10], na.rm = TRUE),
      dRn      = Fluor - baseline
    ) %>%
    summarise(
      Cq = {
        idx <- which(!is.na(dRn) & dRn >= Threshold)[1]
        if (is.na(idx))    NA_real_
        else if (idx == 1) Cycle[1]
        else {
          x0 <- Cycle[idx-1]; y0 <- dRn[idx-1]
          x1 <- Cycle[idx];   y1 <- dRn[idx]
          if (is.na(y0) || is.na(y1) || y1 == y0) NA_real_
          else x0 + (Threshold - y0) / (y1 - y0)
        }
      },
      .groups = "drop"
    ) %>%
    mutate(Threshold = Threshold)
}

# ============================================================================
# SECTION 3: BUILD WELL MAP & RAW Cq VALUES
# ============================================================================

well_map <- pmap_dfr(plate_3slot_design, make_3slot_map_one) %>%
  add_analysis_metadata()

all_cq_raw <- pmap_dfr(plate_files, calc_cq_one_plate)

write_csv(
  all_cq_raw %>% left_join(well_map, by = c("Plate","Well","Row","Col")),
  "all_cq_values_assigned_wells_metadata.csv"
)
cat("Cq values written to CSV\n")

# ============================================================================
# SECTION 4: EXCLUSIONS
# ============================================================================

excluded_wells <- tribble(
  ~Plate,    ~Well,
  "plate13", "G6",
  "plate1",  "G8",
  "plate1",  "G2",
  "plate1",  "G5"
)

exclude_by_tech_rep <- tribble(
  ~Plate,   ~Tech_Rep,
  "plate7", "T3",
  "plate6", "T3"
)

tech_rep_excluded_wells <- all_cq_raw %>%
  left_join(well_map, by = c("Plate","Well","Row","Col")) %>%
  semi_join(exclude_by_tech_rep, by = c("Plate","Tech_Rep")) %>%
  distinct(Plate, Well)

bio_rep_excluded_wells <- all_cq_raw %>%
  left_join(well_map, by = c("Plate","Well","Row","Col")) %>%
  filter(
    (Plate == "plate10" & Condition == "D50"      & Bio_Rep == "R3") |
    (Plate == "plate9"  & Condition == "WT K562"  & Bio_Rep == "R3") |
    (Plate == "plate12" & Condition == "aG ctrl"  & Bio_Rep == "R3")
  ) %>%
  distinct(Plate, Well)

missing_cq_wells <- all_cq_raw %>%
  filter(is.na(Cq)) %>%
  distinct(Plate, Well)

excluded_wells_final <- bind_rows(
  excluded_wells, tech_rep_excluded_wells,
  bio_rep_excluded_wells, missing_cq_wells
) %>% distinct()

# ============================================================================
# SECTION 5: NORMALIZATION (ΔCq → ΔΔCq → fold change)
# ============================================================================

mean_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
sd_na   <- function(x) if (sum(!is.na(x)) <= 1) NA_real_ else sd(x, na.rm = TRUE)
se_na   <- function(x) { n <- sum(!is.na(x)); if (n <= 1) NA_real_ else sd(x, na.rm = TRUE) / sqrt(n) }

get_control_mean <- function(tbl, plate, target_gene, condition) {
  out <- tbl %>%
    filter(Plate == plate, TargetGene == target_gene, Condition == condition) %>%
    pull(control_mean_dCq)
  if (length(out) == 0) NA_real_ else out[1]
}

sample_cols_tech <- c(
  "Plate","Slot","System","CellLine","Guide","TargetGene",
  "Condition","DoseLabel","DoseType","DoseValue",
  "ControlGroup","IsControl","NormalizeToCondition","NormalizeToGroup",
  "Bio_Rep","Tech_Rep"
)

sample_cols_no_tech <- setdiff(sample_cols_tech, "Tech_Rep")

cq_raw <- all_cq_raw %>%
  left_join(well_map, by = c("Plate","Well","Row","Col")) %>%
  filter(!is.na(Sample), !is.na(Gene)) %>%
  anti_join(excluded_wells_final, by = c("Plate","Well"))

# Step 1: paired technical-replicate ΔCq
cq_target_tech <- cq_raw %>%
  filter(Gene == TargetGene) %>%
  select(all_of(sample_cols_tech), target_Cq = Cq)

cq_hk_tech <- cq_raw %>%
  filter(Gene %in% c("Actin","18S")) %>%
  group_by(across(all_of(sample_cols_tech))) %>%
  summarise(
    Actin_Cq           = mean_na(Cq[Gene == "Actin"]),
    `18S_Cq`           = mean_na(Cq[Gene == "18S"]),
    housekeeper_mean_Cq = mean_na(Cq),
    .groups = "drop"
  )

cq_dCq_tech <- cq_target_tech %>%
  left_join(cq_hk_tech, by = sample_cols_tech) %>%
  mutate(dCq_tech = target_Cq - housekeeper_mean_Cq)

# Step 2: average across technical replicates → biological-rep ΔCq
cq_dCq <- cq_dCq_tech %>%
  group_by(across(all_of(sample_cols_no_tech))) %>%
  summarise(
    dCq             = mean_na(dCq_tech),
    target_mean_Cq  = mean_na(target_Cq),
    housekeeper_mean_Cq = mean_na(housekeeper_mean_Cq),
    n_tech          = sum(!is.na(dCq_tech)),
    .groups = "drop"
  )

# Step 3: control means table (used for NTC comparison)
control_means <- cq_dCq %>%
  filter(IsControl) %>%
  group_by(Plate, TargetGene, Condition, ControlGroup, CellLine, Guide) %>%
  summarise(
    control_mean_dCq = mean_na(dCq),
    control_sd_dCq   = sd_na(dCq),
    control_n_bio    = sum(!is.na(dCq)),
    .groups = "drop"
  )

# Step 4: plate-matched control means
plate_matched_controls <- cq_dCq %>%
  filter(IsControl, Condition == NormalizeToCondition) %>%
  group_by(Plate, TargetGene, NormalizeToCondition, NormalizeToGroup) %>%
  summarise(matched_ctrl_mean_dCq = mean_na(dCq), .groups = "drop")

# Step 4: ΔΔCq and fold change
cq_no_bridge <- cq_dCq %>%
  left_join(plate_matched_controls,
            by = c("Plate","TargetGene","NormalizeToCondition","NormalizeToGroup")) %>%
  mutate(
    AnalysisGroup    = paste(CellLine, TargetGene),
    ddCq_no_bridge   = dCq - matched_ctrl_mean_dCq,
    log2fc_no_bridge = -ddCq_no_bridge,
    fold_change_no_bridge = 2^log2fc_no_bridge
  )

# Create plotting objects used by the plot functions
convert_dose_to_uM <- function(DoseType, DoseValue) {
  case_when(
    DoseType == "Dox" ~ as.numeric(DoseValue) / 1000,
    DoseType == "TMP" ~ as.numeric(DoseValue),
    TRUE              ~ as.numeric(DoseValue)
  )
}

no_bridge_plot_data <- cq_no_bridge %>%
  filter(DoseType %in% c("Dox","TMP"), !is.na(DoseValue)) %>%
  mutate(
    concentration = convert_dose_to_uM(DoseType, DoseValue),
    log2fc        = log2fc_no_bridge,
    fold_change   = fold_change_no_bridge
  )

no_bridge_means <- no_bridge_plot_data %>%
  group_by(CellLine, TargetGene, DoseType, DoseValue, concentration) %>%
  summarise(
    mean_log2fc = mean_na(log2fc),
    se_log2fc   = se_na(log2fc),
    mean_fc     = mean_na(fold_change),
    se_fc       = se_na(fold_change),
    n_bio       = sum(!is.na(log2fc)),
    .groups     = "drop"
  )

cat("Normalization complete\n")

# ============================================================================
# SECTION 6: NTC COMPARISON (plate8/9 vs iNTC, iC, and aC controls)
# ============================================================================

plate8_iC_CD46   <- get_control_mean(control_means, "plate8",  "CD46", "iC ctrl")
plate8_iNTC_CD46 <- get_control_mean(control_means, "plate8",  "CD46", "iNTC ctrl")
plate9_iNTC_CD46 <- get_control_mean(control_means, "plate9",  "CD46", "iNTC ctrl")
plate14_iC_CD46  <- get_control_mean(control_means, "plate14", "CD46", "iC ctrl")
plate14_aC_CD46  <- get_control_mean(control_means, "plate14", "CD46", "aC ctrl")

ntc_plate_offsets <- tibble(
  Plate               = c("plate8", "plate9"),
  ntc_anchor_dCq      = c(plate8_iNTC_CD46, plate9_iNTC_CD46),
  plate8_iNTC_ref_dCq = plate8_iNTC_CD46,
  plate8_iC_ref_dCq   = plate8_iC_CD46
) %>%
  mutate(plate_offset_to_plate8_iNTC = ntc_anchor_dCq - plate8_iNTC_ref_dCq)

cq_ntc_comparison <- cq_dCq %>%
  filter(Plate %in% c("plate8", "plate9")) %>%
  left_join(ntc_plate_offsets, by = "Plate") %>%
  mutate(
    dCq_plate8_iNTC_scale = dCq - plate_offset_to_plate8_iNTC,
    log2fc_vs_iNTC        = -(dCq_plate8_iNTC_scale - plate8_iNTC_ref_dCq),
    log2fc_vs_iC          = -(dCq_plate8_iNTC_scale - plate8_iC_ref_dCq),
    dCq_plate14_iC_scale  = dCq_plate8_iNTC_scale - (plate8_iC_CD46 - plate14_iC_CD46),
    log2fc_vs_aC          = -(dCq_plate14_iC_scale - plate14_aC_CD46)
  )

write_csv(cq_ntc_comparison, "plate8_9_ntc_comparison.csv")

p_ntc_comparison <- cq_ntc_comparison %>%
  select(Plate, Condition, Bio_Rep, log2fc_vs_iNTC, log2fc_vs_iC, log2fc_vs_aC) %>%
  pivot_longer(cols = starts_with("log2fc_vs_"), names_to = "Reference", values_to = "log2fc") %>%
  mutate(Reference = recode(Reference,
    "log2fc_vs_iNTC" = "vs iNTC ctrl",
    "log2fc_vs_iC"   = "vs iC ctrl",
    "log2fc_vs_aC"   = "vs aC ctrl"
  )) %>%
  ggplot(aes(Condition, log2fc, color = Condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(shape = Plate),
             position = position_jitter(width = 0.12, height = 0, seed = 42),
             size = 2.8, alpha = 0.75) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, color = "black") +
  facet_wrap(~ Reference, scales = "free_y") +
  labs(x = "Condition", y = bquote(log[2] ~ fold ~ change),
       color = "Condition", shape = "Plate") +
  theme_nathan() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================================
# SECTION 7: PLOT FUNCTIONS
# ============================================================================

plot_cd46_hill <- function(cell_line_choice, dose_type_choice, exclude_fit_conc = NULL) {
  plot_color <- unname(system_colors[cell_line_choice])
  raw_color  <- scales::alpha(plot_color, 0.35)

  plot_data <- no_bridge_plot_data %>%
    filter(CellLine == cell_line_choice, TargetGene == "CD46",
           DoseType == dose_type_choice, !is.na(concentration),
           concentration > 0, !is.na(log2fc))

  if (nrow(plot_data) == 0) stop("No data found.")

  mean_data <- plot_data %>%
    group_by(concentration) %>%
    summarise(
      mean_log2fc = mean_na(log2fc),
      se_log2fc   = se_na(log2fc),
      .groups = "drop"
    ) %>%
    arrange(concentration)

  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.7) +
    geom_point(data = plot_data, aes(concentration, log2fc),
               alpha = 0.45, size = 2.5, colour = raw_color) +
    geom_errorbar(data = mean_data,
                  aes(concentration, ymin = mean_log2fc - se_log2fc, ymax = mean_log2fc + se_log2fc),
                  width = 0.08, linewidth = 0.8, colour = "black") +
    geom_point(data = mean_data, aes(concentration, mean_log2fc),
               size = 4.5, colour = plot_color) +
    scale_x_log10(breaks = sort(unique(mean_data$concentration)),
                  labels = sort(unique(mean_data$concentration))) +
    labs(x = paste0(dose_type_choice, " concentration (\u00b5M, log scale)"),
         y = bquote(CD46 ~ log[2] ~ fold ~ change)) +
    theme_nathan() +
    theme(legend.position = "none")

  # Hill fit
  fit_data <- mean_data
  if (!is.null(exclude_fit_conc)) fit_data <- filter(fit_data, !concentration %in% exclude_fit_conc)
  if (dose_type_choice == "TMP") {
    fit_data <- bind_rows(tibble(concentration = 1e-8, mean_log2fc = 0), fit_data) %>%
      arrange(concentration)
  }

  if (nrow(fit_data) >= 4) {
    low  <- mean_na(fit_data$mean_log2fc[fit_data$concentration == min(fit_data$concentration)])
    high <- mean_na(fit_data$mean_log2fc[fit_data$concentration == max(fit_data$concentration)])
    direction    <- if (high > low) "activation" else "inhibition"
    bottom_fixed <- min(fit_data$mean_log2fc, na.rm = TRUE)
    top_fixed    <- if (direction == "activation") max(fit_data$mean_log2fc, na.rm = TRUE) else 0

    hill_model <- tryCatch(
      drc::drm(mean_log2fc ~ concentration, data = fit_data,
               fct = drc::LL.4(fixed = c(NA, bottom_fixed, top_fixed, NA),
                               names = c("Hill","Bottom","Top","EC50"))),
      error = function(e) { message("Hill fit failed: ", e$message); NULL }
    )

    if (!is.null(hill_model)) {
      min_c <- if (dose_type_choice == "TMP") 1e-8 else min(fit_data$concentration[fit_data$concentration > 0])
      conc_range <- data.frame(concentration = exp(seq(log(min_c), log(max(fit_data$concentration)), length.out = 500)))

      pred <- suppressWarnings(tryCatch(
        predict(hill_model, newdata = conc_range, interval = "confidence", level = 0.95),
        error = function(e) NULL
      ))

      pred_hill <- if (!is.null(pred)) {
        conc_range %>% mutate(fit = pred[,1], ci_lo = pred[,2], ci_hi = pred[,3])
      } else {
        conc_range %>% mutate(
          fit   = suppressWarnings(tryCatch(as.numeric(predict(hill_model, newdata = conc_range)), error = function(e) NA_real_)),
          ci_lo = NA_real_, ci_hi = NA_real_
        )
      }
      pred_hill <- filter(pred_hill, is.finite(fit))

      if (nrow(pred_hill) > 0) {
        if (!all(is.na(pred_hill$ci_lo))) {
          p <- p + geom_ribbon(data = pred_hill, aes(concentration, ymin = ci_lo, ymax = ci_hi),
                               fill = plot_color, alpha = 0.18, inherit.aes = FALSE)
        }
        p <- p + geom_line(data = pred_hill, aes(concentration, fit),
                           linewidth = 1.1, colour = plot_color)
      }
    }
  } else {
    message("Skipping Hill fit — fewer than 4 dose points.")
  }
  p
}

plot_gfi1b_bar <- function(cell_line_choice, dose_type_choice) {
  plot_color <- unname(system_colors[cell_line_choice])

  plot_data <- no_bridge_plot_data %>%
    filter(CellLine == cell_line_choice, TargetGene == "GFI1B",
           DoseType == dose_type_choice, !is.na(concentration), concentration > 0)

  if (nrow(plot_data) == 0) stop("No GFI1B data found.")

  mean_data <- plot_data %>%
    group_by(concentration) %>%
    summarise(
      n      = sum(!is.na(log2fc)),
      mean_y = mean_na(log2fc),
      ci_y   = qt(0.975, df = pmax(n-1, 1)) * se_na(log2fc),
      .groups = "drop"
    ) %>%
    arrange(concentration)

  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.7) +
    geom_col(data = mean_data, aes(factor(concentration), mean_y),
             fill = plot_color, alpha = 0.65, width = 0.6) +
    geom_errorbar(data = mean_data,
                  aes(factor(concentration), ymin = mean_y - ci_y, ymax = mean_y + ci_y),
                  width = 0.18, linewidth = 0.8, color = "black") +
    geom_point(data = plot_data,
               aes(factor(concentration), log2fc),
               color = plot_color,
               position = position_jitter(width = 0.08, height = 0, seed = 42),
               size = 3, alpha = 0.85) +
    labs(x = paste0(dose_type_choice, " concentration (\u00b5M)"),
         y = bquote(GFI1B ~ log[2] ~ fold ~ change)) +
    theme_nathan() +
    theme(legend.position = "none")
}

# ============================================================================
# SECTION 8: GENERATE & SAVE ALL PLOTS
# ============================================================================

p_i_cd46_dox  <- plot_cd46_hill("CRISPRi", "Dox")
p_i_cd46_tmp  <- plot_cd46_hill("CRISPRi", "TMP", exclude_fit_conc = 10)
p_a_cd46_dox  <- plot_cd46_hill("CRISPRa", "Dox")

p_bar_i_gfi1b_dox <- plot_gfi1b_bar("CRISPRi", "Dox")
p_bar_i_gfi1b_tmp <- plot_gfi1b_bar("CRISPRi", "TMP")
p_bar_a_gfi1b_dox <- plot_gfi1b_bar("CRISPRa", "Dox")

target_plots <- list(
  p_i_cd46_dox      = p_i_cd46_dox,
  p_i_cd46_tmp      = p_i_cd46_tmp,
  p_a_cd46_dox      = p_a_cd46_dox,
  p_bar_i_gfi1b_dox = p_bar_i_gfi1b_dox,
  p_bar_i_gfi1b_tmp = p_bar_i_gfi1b_tmp,
  p_bar_a_gfi1b_dox = p_bar_a_gfi1b_dox
)

walk2(names(target_plots), target_plots, function(name, plot) {
  ggsave(file.path(output_dir, paste0(name, ".png")), plot, width = 5, height = 4, dpi = 300)
  ggsave(file.path(output_dir, paste0(name, ".pdf")), plot, width = 5, height = 4)
  cat("Saved:", name, "\n")
})

ggsave(file.path(output_dir, "p_ntc_comparison.png"), p_ntc_comparison, width = 9, height = 5, dpi = 300)
ggsave(file.path(output_dir, "p_ntc_comparison.pdf"), p_ntc_comparison, width = 9, height = 5)
cat("Saved: p_ntc_comparison\n")
