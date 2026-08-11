#!/usr/bin/env Rscript

# Read-only statistical audit for the current Figure 3 and Figure 4 outputs.
# Produces alternative multiplicity columns and fit-QC summaries without
# changing any figure, threshold, test, or manuscript statement.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out_dir <- file.path(project_root, "results")

fig3_file <- file.path(out_dir, "Figure3_source_data_v16.xlsx")
fig4_file <- file.path(out_dir, "Fig4_Data_V6.xlsx")
stopifnot(file.exists(fig3_file), file.exists(fig4_file))

add_adjustments <- function(data, endpoint, p_col = "p_value") {
  data %>%
    mutate(
      endpoint = endpoint,
      p_nominal = .data[[p_col]],
      p_holm = p.adjust(p_nominal, method = "holm"),
      p_bh = p.adjust(p_nominal, method = "BH"),
      sig_nominal = is.finite(p_nominal) & p_nominal < 0.05,
      sig_holm = is.finite(p_holm) & p_holm < 0.05,
      sig_bh = is.finite(p_bh) & p_bh < 0.05
    )
}

fig3_audit <- bind_rows(
  add_adjustments(read_excel(fig3_file, sheet = "ttest_t50_per_dose"), "t50"),
  add_adjustments(read_excel(fig3_file, sheet = "ttest_rel_per_dose"), "relative germination fraction at 3 h")
) %>%
  relocate(endpoint, Dose, p_nominal, p_holm, p_bh,
           sig_nominal, sig_holm, sig_bh)

write_csv(fig3_audit, file.path(out_dir, "Figure3_multiplicity_audit_v18.csv"))

fit_metrics <- read_excel(fig4_file, sheet = "fitted_metrics")
fig4_fit_qc <- fit_metrics %>%
  mutate(
    fit_ok = as.logical(fit_ok),
    hit_bound = as.logical(hit_bound),
    qc_status = case_when(
      !fit_ok ~ "fit failed",
      hit_bound ~ "parameter at bound",
      !is.finite(mu_max) | !is.finite(lag_time) ~ "non-finite derived metric",
      TRUE ~ "passed current checks"
    )
  )

write_csv(fig4_fit_qc, file.path(out_dir, "Figure4_fit_QC_rows_v18.csv"))

fig4_fit_qc_summary <- fig4_fit_qc %>%
  count(Condition, Dose, qc_status, name = "n") %>%
  complete(Condition, Dose, qc_status, fill = list(n = 0))

write_csv(fig4_fit_qc_summary, file.path(out_dir, "Figure4_fit_QC_summary_v18.csv"))

fig4_emm_audit <- bind_rows(
  add_adjustments(read_excel(fig4_file, sheet = "emmeans_mu"), "mu_max", "p.value"),
  add_adjustments(read_excel(fig4_file, sheet = "emmeans_lag"), "lag_time", "p.value")
) %>%
  relocate(endpoint, Dose_f, contrast, p_nominal, p_holm, p_bh,
           sig_nominal, sig_holm, sig_bh)

write_csv(fig4_emm_audit, file.path(out_dir, "Figure4_multiplicity_audit_v18.csv"))

message("Figure 3/4 statistical audit completed without changing submitted figures.")
