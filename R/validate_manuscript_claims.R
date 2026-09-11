project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out_dir <- file.path(project_root, "results")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

one_value <- function(data, column) {
  if (nrow(data) != 1L) {
    stop("Expected exactly one matching row, found ", nrow(data), call. = FALSE)
  }
  as.numeric(data[[column]][[1]])
}

checks <- tibble(
  claim = character(),
  observed = numeric(),
  expected = numeric(),
  tolerance = numeric()
)

add_check <- function(claim, observed, expected, tolerance) {
  checks <<- bind_rows(
    checks,
    tibble(
      claim = claim,
      observed = as.numeric(observed),
      expected = as.numeric(expected),
      tolerance = as.numeric(tolerance)
    )
  )
}

# Inoculum-uniformity statements in Methods and Supplementary Information.
inoculum <- read_csv(file.path(out_dir, "inoculum_OD600_uniformity_qc.csv"), show_col_types = FALSE)
add_check(
  "Overall pre-transfer OD600",
  one_value(filter(inoculum, metric == "overall_mean_OD600_0h"), "value"),
  0.0357034717083333,
  1e-10
)
add_check(
  "CV among eight treatment-group OD600 means (%)",
  one_value(filter(inoculum, metric == "CV_percent_among_eight_treatment_group_means"), "value"),
  4.30413403667593,
  1e-8
)

# Stage-resolved germination claims.
fig3 <- read_csv(file.path(out_dir, "Figure3_multiplicity_audit_v18.csv"), show_col_types = FALSE)
for (spec in list(
  list("t50: SMG minus 1 g at 6.5 Gy (h)", "t50", 6.5, "estimate_SMG_minus_1g", -0.418223789080057),
  list("t50: Holm p at 6.5 Gy", "t50", 6.5, "p_holm", 0.0439217686547408),
  list("Relative germination: SMG minus 1 g at 13 Gy", "relative germination fraction at 3 h", 13, "estimate_SMG_minus_1g", -0.209095280828411),
  list("Relative germination: Holm p at 13 Gy", "relative germination fraction at 3 h", 13, "p_holm", 0.003403852443036752),
  list("Relative germination: SMG minus 1 g at 19.5 Gy", "relative germination fraction at 3 h", 19.5, "estimate_SMG_minus_1g", -0.146510301749208),
  list("Relative germination: Holm p at 19.5 Gy", "relative germination fraction at 3 h", 19.5, "p_holm", 0.030559814774007602)
)) {
  row <- filter(fig3, endpoint == spec[[2]], Dose == spec[[3]])
  add_check(spec[[1]], one_value(row, spec[[4]]), spec[[5]], 1e-8)
}

# Post-germination growth claim. The audit reports 1 g minus SMG.
fig4 <- read_csv(file.path(out_dir, "Figure4_multiplicity_audit_v18.csv"), show_col_types = FALSE)
r_65 <- filter(
  fig4,
  endpoint == "apparent population-level growth-rate coefficient r",
  Dose_f == 6.5
)
add_check("Growth-rate coefficient: 1 g minus SMG at 6.5 Gy (h^-1)", one_value(r_65, "estimate"), 0.522383162404782, 1e-8)
add_check("Growth-rate coefficient: Holm p at 6.5 Gy", one_value(r_65, "p_holm"), 0.00701609819568196, 1e-8)

# Transcriptome-level variance and dispersion claims.
permanova <- read_csv(file.path(out_dir, "Figure5_PERMANOVA_reproduced.csv"), show_col_types = FALSE)
for (spec in list(
  list("PERMANOVA 3 h dose R2", "3 h", "dose_f", "R2", 0.20291021356361585),
  list("PERMANOVA 3 h dose p", "3 h", "dose_f", "Pr(>F)", 0.017),
  list("PERMANOVA 3 h gravity R2", "3 h", "condition", "R2", 0.08239043375990066),
  list("PERMANOVA 3 h gravity p", "3 h", "condition", "Pr(>F)", 0.035),
  list("PERMANOVA 10 h dose R2", "10 h", "dose_f", "R2", 0.2217361742023482),
  list("PERMANOVA 10 h dose p", "10 h", "dose_f", "Pr(>F)", 0.001),
  list("PERMANOVA 10 h gravity R2", "10 h", "condition", "R2", 0.09232136712964457),
  list("PERMANOVA 10 h gravity p", "10 h", "condition", "Pr(>F)", 0.006),
  list("PERMANOVA 10 h interaction R2", "10 h", "condition:dose_f", "R2", 0.1762995424724566),
  list("PERMANOVA 10 h interaction p", "10 h", "condition:dose_f", "Pr(>F)", 0.017)
)) {
  row <- filter(permanova, timepoint == spec[[2]], term == spec[[3]])
  add_check(spec[[1]], one_value(row, spec[[4]]), spec[[5]], 1e-8)
}

permdisp <- read_csv(file.path(out_dir, "Figure5_PERMDISP_audit.csv"), show_col_types = FALSE)
add_check("PERMDISP 3 h p", one_value(filter(permdisp, timepoint == "3 h"), "p_value"), 0.919, 1e-8)
add_check("PERMDISP 10 h p", one_value(filter(permdisp, timepoint == "10 h"), "p_value"), 0.905, 1e-8)

# Functional-category/phenotype associations reported in Results.
fig6 <- read_csv(file.path(out_dir, "Figure6_module_phenotype_correlations.csv"), show_col_types = FALSE)
for (spec in list(
  list("Transport versus delta r: Spearman rho", "delta_growth_rate_coef", "Transport", "cor_val", -0.7714285714285715),
  list("Transport versus delta r: n", "delta_growth_rate_coef", "Transport", "n_valid", 6),
  list("Transport versus delta AUC: Spearman rho", "delta_AUC", "Transport", "cor_val", -0.7142857142857143),
  list("Metabolism versus delta r: Spearman rho", "delta_growth_rate_coef", "Metabolism", "cor_val", -0.6071428571428571),
  list("Metabolism versus delta AUC: Spearman rho", "delta_AUC", "Metabolism", "cor_val", -0.6071428571428571)
)) {
  row <- filter(fig6, Phenotype == spec[[2]], Bin == spec[[3]])
  add_check(spec[[1]], one_value(row, spec[[4]]), spec[[5]], 1e-8)
}

# Secondary-spore AUC statements.
fig7 <- read_csv(file.path(out_dir, "Figure7_pooled_baseline_descriptive_v18.csv"), show_col_types = FALSE)
baseline_row <- function(group_name) filter(fig7, .data$run == "H2O2", .data$group == group_name)
challenge_row <- function(run_name, group_name) filter(fig7, .data$run == run_name, .data$group == group_name)
add_check("Secondary spores: SMG + 13 Gy baseline AUC", one_value(baseline_row("SMG+13"), "baseline_mean"), 0.182, 5e-4)
add_check("Secondary spores: SMG + 19.5 Gy baseline AUC", one_value(baseline_row("SMG+19.5"), "baseline_mean"), -0.096, 5e-4)
add_check("Secondary spores: SMG + 13 Gy H2O2 net AUC", one_value(challenge_row("H2O2", "SMG+13"), "descriptive_net_kill"), 0.553, 5e-4)
add_check("Secondary spores: control H2O2 net AUC", one_value(challenge_row("H2O2", "Ctrl"), "descriptive_net_kill"), 0.436, 5e-4)
add_check("Secondary spores: SMG + 13 Gy gentamicin net AUC", one_value(challenge_row("Gentamicin", "SMG+13"), "descriptive_net_kill"), 0.395, 5e-4)
add_check("Secondary spores: 1 g + 13 Gy gentamicin net AUC", one_value(challenge_row("Gentamicin", "1g+13"), "descriptive_net_kill"), 0.072, 5e-4)

# Supplementary-table completeness and curve-fit claims.
table_s1 <- read_csv(file.path(out_dir, "Table_S1_DEG_list.csv"), show_col_types = FALSE)
table_s2 <- read_csv(file.path(out_dir, "Table_S2_module_gene_membership.csv"), show_col_types = FALSE)
table_s3 <- read_csv(file.path(out_dir, "Table_S3_logistic_model_fit_quality.csv"), show_col_types = FALSE)
s3_source <- read_csv(file.path(out_dir, "Supplementary_Figure_S3_source_data.csv"), show_col_types = FALSE)
add_check("Table S1 data rows", nrow(table_s1), 52413, 0)
add_check("Table S1 timepoint-by-contrast combinations", n_distinct(interaction(table_s1$timepoint, table_s1$contrast)), 14, 0)
add_check("Table S2 unique genes", n_distinct(table_s2$gene_id), 181, 0)
add_check("Table S3 replicate fits", nrow(table_s3), 24, 0)
add_check("Table S3 converged fits", sum(table_s3$converged), 24, 0)
add_check("Table S3 mean R2", mean(table_s3$R2), 0.981, 5e-4)
add_check("Table S3 minimum R2", min(table_s3$R2), 0.961, 5e-4)
add_check("Table S3 maximum R2", max(table_s3$R2), 0.993, 5e-4)
add_check("Figure S3 source-data rows", nrow(s3_source), 5040, 0)
add_check("Figure S3 observed points", sum(s3_source$series == "observed"), 216, 0)
add_check("Figure S3 fitted coordinates", sum(s3_source$series == "fitted"), 4824, 0)
add_check("Figure S3 source-data samples", n_distinct(s3_source$SampleID), 24, 0)
add_check("Figure S3 missing OD600 values", sum(is.na(s3_source$OD600)), 0, 0)

checks <- checks %>%
  mutate(
    absolute_difference = abs(observed - expected),
    pass = is.finite(observed) & absolute_difference <= tolerance
  )

write_csv(checks, file.path(out_dir, "manuscript_claim_validation.csv"))

if (any(!checks$pass)) {
  print(filter(checks, !pass), n = Inf)
  stop("Manuscript-claim validation failed; see results/manuscript_claim_validation.csv", call. = FALSE)
}

message("Validated ", nrow(checks), " manuscript and Supplementary Information claims.")
