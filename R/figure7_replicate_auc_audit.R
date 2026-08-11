#!/usr/bin/env Rscript

# Figure 7 replicate-level AUC and baseline-pairing audit.
#
# This script is diagnostic. It does not overwrite the submitted Figure 7 or
# assume that control_1/control_2 are paired to a particular biocide run.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(stringr)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
source_dir <- file.path(project_root, "data", "secondary_spores")
out_dir <- file.path(project_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- c(
  ctrl1 = "control_1.xlsx",
  ctrl2 = "control_2.xlsx",
  H2O2 = "h2o2.xlsx",
  Gentamicin = "gentamicin.xlsx",
  DDAC = "DDAC.xlsx",
  NaHCO3 = "soda_1.xlsx"
)
stopifnot(all(file.exists(file.path(source_dir, files))))

group_map <- c(
  "Ctrl" = "Ctrl",
  "SMG" = "SMG",
  "SR（6.5 Gy）" = "1g+6.5",
  "SR（13 Gy）" = "1g+13",
  "SR（19.5 Gy）" = "1g+19.5",
  "SMG+SR（6.5 Gy）" = "SMG+6.5",
  "SMG+SR（13 Gy）" = "SMG+13",
  "SMG+SR（19.5 Gy）" = "SMG+19.5"
)

trapz_normalized <- function(time, value) {
  keep <- is.finite(time) & is.finite(value)
  time <- time[keep]
  value <- value[keep]
  ord <- order(time)
  time <- time[ord]
  value <- value[ord]
  if (length(time) < 2L || max(time) == min(time)) return(NA_real_)
  sum(diff(time) * (head(value, -1) + tail(value, -1)) / 2) /
    (max(time) - min(time))
}

read_run <- function(run_id, filename) {
  raw <- read_excel(file.path(source_dir, filename)) %>%
    mutate(across(-time, as.numeric))

  raw %>%
    pivot_longer(-time, names_to = "sample", values_to = "od600") %>%
    mutate(
      replicate = as.integer(str_extract(sample, "[0-9]+$")),
      group_raw = str_remove(sample, "-[0-9]+$"),
      group = unname(group_map[group_raw]),
      run = run_id
    ) %>%
    filter(!is.na(group), !is.na(replicate)) %>%
    group_by(run, group, replicate) %>%
    arrange(time, .by_group = TRUE) %>%
    mutate(
      od0 = first(od600),
      relative_od_change = (od600 - od0) / od0
    ) %>%
    ungroup()
}

all_long <- bind_rows(Map(read_run, names(files), unname(files)))

auc_rep <- all_long %>%
  group_by(run, group, replicate) %>%
  summarise(
    n_timepoints = sum(is.finite(relative_od_change)),
    auc = trapz_normalized(time, relative_od_change),
    .groups = "drop"
  )

write_csv(
  auc_rep,
  file.path(out_dir, "Figure7_replicate_AUC_all_v18.csv")
)

auc_summary <- auc_rep %>%
  group_by(run, group) %>%
  summarise(
    n = sum(is.finite(auc)),
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE),
    se_auc = sd_auc / sqrt(n),
    .groups = "drop"
  )

write_csv(
  auc_summary,
  file.path(out_dir, "Figure7_replicate_AUC_summary_v18.csv")
)

# Candidate paired differences are shown for both possible baseline runs. They
# are a sensitivity analysis, not a declaration of the experimental pairing.
biocide_runs <- c("H2O2", "Gentamicin", "DDAC", "NaHCO3")
candidate_pairs <- tidyr::crossing(
  baseline_run = c("ctrl1", "ctrl2"),
  biocide_run = biocide_runs
)

paired_sensitivity <- bind_rows(lapply(seq_len(nrow(candidate_pairs)), function(i) {
  b <- candidate_pairs$baseline_run[i]
  t <- candidate_pairs$biocide_run[i]
  auc_rep %>%
    filter(run == b) %>%
    select(group, replicate, auc_baseline = auc) %>%
    inner_join(
      auc_rep %>%
        filter(run == t) %>%
        select(group, replicate, auc_biocide = auc),
      by = c("group", "replicate")
    ) %>%
    mutate(
      baseline_candidate = b,
      biocide = t,
      net_kill = auc_baseline - auc_biocide
    )
}))

write_csv(
  paired_sensitivity,
  file.path(out_dir, "Figure7_baseline_pairing_sensitivity_replicates_v18.csv")
)

pairing_summary <- paired_sensitivity %>%
  group_by(biocide, baseline_candidate, group) %>%
  summarise(
    n_pairs = sum(is.finite(net_kill)),
    mean_net_kill = mean(net_kill, na.rm = TRUE),
    sd_net_kill = sd(net_kill, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(biocide, baseline_candidate, group, n_pairs, mean_net_kill, sd_net_kill) %>%
  pivot_wider(
    names_from = baseline_candidate,
    values_from = c(n_pairs, mean_net_kill, sd_net_kill)
  ) %>%
  mutate(
    sign_agrees = case_when(
      !is.finite(mean_net_kill_ctrl1) | !is.finite(mean_net_kill_ctrl2) ~ NA,
      mean_net_kill_ctrl1 == 0 | mean_net_kill_ctrl2 == 0 ~ TRUE,
      TRUE ~ sign(mean_net_kill_ctrl1) == sign(mean_net_kill_ctrl2)
    ),
    pairing_range = abs(mean_net_kill_ctrl1 - mean_net_kill_ctrl2)
  )

write_csv(
  pairing_summary,
  file.path(out_dir, "Figure7_baseline_pairing_sensitivity_summary_v18.csv")
)

# Pooled-baseline contrast is included only as an unpaired descriptive option.
pooled_baseline <- auc_rep %>%
  filter(run %in% c("ctrl1", "ctrl2")) %>%
  group_by(group) %>%
  summarise(
    baseline_n = sum(is.finite(auc)),
    baseline_mean = mean(auc, na.rm = TRUE),
    baseline_sd = sd(auc, na.rm = TRUE),
    .groups = "drop"
  )

pooled_contrast <- auc_rep %>%
  filter(run %in% biocide_runs) %>%
  group_by(run, group) %>%
  summarise(
    biocide_n = sum(is.finite(auc)),
    biocide_mean = mean(auc, na.rm = TRUE),
    biocide_sd = sd(auc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(pooled_baseline, by = "group") %>%
  mutate(descriptive_net_kill = baseline_mean - biocide_mean)

write_csv(
  pooled_contrast,
  file.path(out_dir, "Figure7_pooled_baseline_descriptive_v18.csv")
)

message("Figure 7 replicate-AUC audit completed without changing the submitted figure.")
