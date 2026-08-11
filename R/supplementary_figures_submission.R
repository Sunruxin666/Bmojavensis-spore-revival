## Reproducible Supplementary Figures S1-S3 and Supplementary Tables S1/S4.
##
## Figure contract
## S1: all eight parental pretreatment groups are shown for the two secondary
##     antimicrobial assays that are not part of the main Figure 7.
## S2: PCA represents all 48 RNA-seq libraries (24 per time point), rather than
##     a subset analysis.
## S3: individual growth-model fits expose replicate-level fit quality behind
##     the summary parameters in Figure 4.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(minpack.lm)
  library(patchwork)
  library(readr)
  library(readxl)
  library(tidyr)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
data_root <- file.path(project_root, "data")
out_dir <- file.path(project_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Generic sans resolves to the platform's Arial/Helvetica-compatible family
# while avoiding PostScript metric warnings on headless validation systems.
font_family <- "sans"
fig_width_cm <- 18.3
dpi <- 600

theme_si <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6.5),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      legend.key.size = grid::unit(3.2, "mm"),
      strip.text = element_text(size = 6.5, face = "bold"),
      strip.background = element_rect(fill = "grey96", colour = "grey75", linewidth = 0.3),
      panel.grid = element_blank(),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.margin = margin(4, 5, 4, 5, unit = "pt")
    )
}

save_bundle <- function(plot, stem, height_cm) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot,
         width = fig_width_cm, height = height_cm, units = "cm", device = cairo_pdf)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(file.path(out_dir, paste0(stem, ".svg")), plot,
           width = fig_width_cm, height = height_cm, units = "cm", device = svglite::svglite)
  }
  ggsave(file.path(out_dir, paste0(stem, ".tiff")), plot,
         width = fig_width_cm, height = height_cm, units = "cm", dpi = dpi,
         device = ragg::agg_tiff, compression = "lzw")
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot,
         width = fig_width_cm, height = height_cm, units = "cm", dpi = dpi,
         device = ragg::agg_png)
}

## -------------------------------------------------------------------------
## Supplementary Figure S1: DDAC and NaHCO3 responses across all eight groups
## -------------------------------------------------------------------------

spore_dir <- file.path(data_root, "secondary_spores")
group_cols <- list(
  "1g+0" = c("Ctrl-1", "Ctrl-2", "Ctrl-3"),
  "1g+6.5" = c("SR（6.5 Gy）-1", "SR（6.5 Gy）-2", "SR（6.5 Gy）-3"),
  "1g+13" = c("SR（13 Gy）-1", "SR（13 Gy）-2", "SR（13 Gy）-3"),
  "1g+19.5" = c("SR（19.5 Gy）-1", "SR（19.5 Gy）-2", "SR（19.5 Gy）-3"),
  "SMG+0" = c("SMG-1", "SMG-2", "SMG-3"),
  "SMG+6.5" = c("SMG+SR（6.5 Gy）-1", "SMG+SR（6.5 Gy）-2", "SMG+SR（6.5 Gy）-3"),
  "SMG+13" = c("SMG+SR（13 Gy）-1", "SMG+SR（13 Gy）-2", "SMG+SR（13 Gy）-3"),
  "SMG+19.5" = c("SMG+SR（19.5 Gy）-1", "SMG+SR（19.5 Gy）-2", "SMG+SR（19.5 Gy）-3")
)
group_levels <- names(group_cols)
group_labels <- c(
  "1g+0" = "1g + 0 Gy", "1g+6.5" = "1g + 6.5 Gy",
  "1g+13" = "1g + 13 Gy", "1g+19.5" = "1g + 19.5 Gy",
  "SMG+0" = "SMG + 0 Gy", "SMG+6.5" = "SMG + 6.5 Gy",
  "SMG+13" = "SMG + 13 Gy", "SMG+19.5" = "SMG + 19.5 Gy"
)
group_colors <- c(
  "1g+0" = "#1F78B4", "1g+6.5" = "#6BAED6", "1g+13" = "#4292C6", "1g+19.5" = "#08519C",
  "SMG+0" = "#F28E2B", "SMG+6.5" = "#F4A261", "SMG+13" = "#D95F02", "SMG+19.5" = "#9C3D10"
)

read_od <- function(path) read_excel(path) |> mutate(across(-time, as.numeric))

wide_to_relative <- function(df_wide, treatment) {
  time_vec <- df_wide$time
  bind_rows(lapply(names(group_cols), function(grp) {
    cols <- intersect(group_cols[[grp]], colnames(df_wide))
    if (length(cols) != 3L) stop("Expected three challenge replicates for ", grp, " in ", treatment)
    sub <- as.data.frame(df_wide[, cols, drop = FALSE])
    od0 <- as.numeric(sub[1, ])
    rel <- sweep(as.matrix(sub), 2, od0, "-")
    rel <- sweep(rel, 2, od0, "/")
    as.data.frame(rel) |>
      mutate(time_min = time_vec) |>
      pivot_longer(-time_min, names_to = "rep", values_to = "relative_OD") |>
      mutate(group = grp, treatment = treatment)
  })) |>
    mutate(group = factor(group, levels = group_levels), time_h = time_min / 60)
}

baseline_long <- bind_rows(
  wide_to_relative(read_od(file.path(spore_dir, "control_1.xlsx")), "baseline_1"),
  wide_to_relative(read_od(file.path(spore_dir, "control_2.xlsx")), "baseline_2")
)
challenge_long <- bind_rows(
  wide_to_relative(read_od(file.path(spore_dir, "DDAC.xlsx")), "DDAC"),
  wide_to_relative(read_od(file.path(spore_dir, "soda_1.xlsx")), "NaHCO3")
)

baseline_auc <- baseline_long |>
  group_by(group, time_min) |>
  summarise(relative_OD = mean(relative_OD, na.rm = TRUE), .groups = "drop") |>
  group_by(group) |>
  arrange(time_min, .by_group = TRUE) |>
  summarise(
    baseline_AUC = sum(diff(time_min) * (head(relative_OD, -1) + tail(relative_OD, -1)) / 2) /
      (max(time_min) - min(time_min)),
    .groups = "drop"
  )

challenge_auc <- challenge_long |>
  group_by(treatment, group, rep) |>
  arrange(time_min, .by_group = TRUE) |>
  summarise(
    challenge_AUC = sum(diff(time_min) * (head(relative_OD, -1) + tail(relative_OD, -1)) / 2) /
      (max(time_min) - min(time_min)),
    .groups = "drop"
  ) |>
  left_join(baseline_auc, by = "group") |>
  mutate(net_AUC = baseline_AUC - challenge_AUC)

kinetics_summary <- challenge_long |>
  group_by(treatment, group, time_h) |>
  summarise(mean = mean(relative_OD), sd = sd(relative_OD), .groups = "drop")

write_csv(challenge_long, file.path(out_dir, "Supplementary_Figure_S1_relative_OD_source_data.csv"))
write_csv(challenge_auc, file.path(out_dir, "Supplementary_Figure_S1_net_AUC_source_data.csv"))

plot_kinetics <- function(treatment_label) {
  dd <- filter(kinetics_summary, treatment == treatment_label)
  display_label <- ifelse(treatment_label == "NaHCO3", "NaHCO₃", treatment_label)
  ggplot(dd, aes(time_h, mean, colour = group, fill = group, group = group)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.10, colour = NA) +
    geom_line(linewidth = 0.55) +
    scale_colour_manual(values = group_colors, labels = group_labels, drop = FALSE) +
    scale_fill_manual(values = group_colors, labels = group_labels, drop = FALSE) +
    labs(x = "Time (h)", y = "Relative OD₆₀₀", title = display_label) +
    theme_si() + theme(legend.position = "none")
}

plot_auc <- function(treatment_label) {
  dd <- filter(challenge_auc, treatment == treatment_label)
  display_label <- ifelse(treatment_label == "NaHCO3", "NaHCO₃", treatment_label)
  ggplot(dd, aes(group, net_AUC, colour = group)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed", colour = "grey55") +
    geom_point(position = position_jitter(width = 0.10, height = 0, seed = 20260718), size = 1.3) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2.0, colour = "black") +
    scale_colour_manual(values = group_colors, labels = group_labels, drop = FALSE) +
    scale_x_discrete(labels = group_labels) +
    labs(x = NULL, y = "Post-challenge net AUC", title = display_label) +
    theme_si() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

legend_plot <- ggplot(
  data.frame(group = factor(group_levels, levels = group_levels), x = seq_along(group_levels), y = 1),
  aes(x, y, colour = group)
) +
  geom_point(size = 2) +
  scale_colour_manual(values = group_colors, labels = group_labels, drop = FALSE) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, title = NULL)) +
  theme_void(base_family = font_family) +
  theme(legend.position = "bottom", legend.text = element_text(size = 6))

shared_legend <- cowplot::get_legend(legend_plot)
s1_main <- (((plot_kinetics("DDAC") + labs(tag = "a")) |
             (plot_kinetics("NaHCO3") + labs(tag = "b"))) /
            ((plot_auc("DDAC") + labs(tag = "c")) |
             (plot_auc("NaHCO3") + labs(tag = "d")))) +
  plot_layout(heights = c(1.05, 1))
s1 <- wrap_plots(
  list(s1_main, wrap_elements(full = shared_legend)),
  ncol = 1,
  heights = c(2.05, 0.18)
)
save_bundle(s1, "Supplementary_Figure_S1", 15.5)

## -------------------------------------------------------------------------
## Supplementary Figure S2: PCA of all RNA-seq samples
## -------------------------------------------------------------------------

count_root <- file.path(data_root, "transcriptomics", "count_matrices")

merge_count_files <- function(timepoint) {
  paths <- file.path(
    count_root, timepoint,
    c("Ctrl_all_CountMatrix.tsv", "Ctrl_vs_SMG_CountMatrix.tsv", "Ctrl_vs_SMG_SR_CountMatrix.tsv")
  )
  if (!all(file.exists(paths))) stop("Missing count matrix for ", timepoint)
  mats <- lapply(paths, read_tsv, show_col_types = FALSE)
  genes <- Reduce(intersect, lapply(mats, function(x) x$Geneid))
  mats <- lapply(mats, function(x) as.data.frame(x[match(genes, x$Geneid), ]))
  out <- data.frame(Geneid = genes, check.names = FALSE)
  sample_names <- unique(unlist(lapply(mats, function(x) setdiff(names(x), "Geneid"))))
  for (sample_name in sample_names) {
    candidates <- lapply(mats, function(x) if (sample_name %in% names(x)) x[[sample_name]] else NULL)
    candidates <- candidates[!vapply(candidates, is.null, logical(1))]
    if (length(candidates) > 1L) {
      identical_values <- all(vapply(candidates[-1], identical, logical(1), candidates[[1]]))
      if (!identical_values) stop("Conflicting duplicated count column: ", sample_name, " at ", timepoint)
    }
    out[[sample_name]] <- as.numeric(candidates[[1]])
  }
  out
}

parse_sample <- function(sample_names, timepoint) {
  tibble(sample = sample_names) |>
    mutate(
      gravity = if_else(grepl("^SMG", sample), "SMG", "1g"),
      dose = case_when(
        grepl("19_5", sample) ~ 19.5,
        grepl("6_5", sample) ~ 6.5,
        grepl("13", sample) ~ 13,
        TRUE ~ 0
      ),
      timepoint = timepoint,
      gravity = factor(gravity, levels = c("1g", "SMG")),
      dose_f = factor(dose, levels = c(0, 6.5, 13, 19.5))
    )
}

run_pca <- function(timepoint) {
  counts <- merge_count_files(timepoint)
  meta <- parse_sample(setdiff(names(counts), "Geneid"), timepoint)
  if (nrow(meta) != 24L) stop("Expected 24 samples at ", timepoint, "; found ", nrow(meta))
  mat <- as.matrix(counts[, meta$sample, drop = FALSE])
  storage.mode(mat) <- "numeric"
  keep <- rowSums(mat) >= 20
  mat <- mat[keep, , drop = FALSE]
  lib_sizes <- colSums(mat)
  log_cpm <- log2(sweep(mat, 2, lib_sizes / 1e6, "/") + 1)
  pca <- prcomp(t(log_cpm), center = TRUE, scale. = FALSE)
  pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  coords <- as.data.frame(pca$x[, 1:2, drop = FALSE]) |>
    tibble::rownames_to_column("sample") |>
    left_join(meta, by = "sample") |>
    mutate(PC1_percent = pct[1], PC2_percent = pct[2])
  qc <- tibble(
    sample = meta$sample,
    timepoint = timepoint,
    gravity = as.character(meta$gravity),
    dose_Gy = meta$dose,
    total_assigned_counts = colSums(mat),
    detected_genes = colSums(mat > 0),
    median_nonzero_count = apply(mat, 2, function(x) median(x[x > 0])),
    library_size_factor = colSums(mat) / exp(mean(log(colSums(mat))))
  )
  list(coords = coords, qc = qc, pct = pct)
}

pca3 <- run_pca("3h")
pca10 <- run_pca("10h")
pca_coords <- bind_rows(pca3$coords, pca10$coords)
qc_table <- bind_rows(pca3$qc, pca10$qc)
stopifnot(nrow(pca_coords) == 48L, nrow(qc_table) == 48L)
write_csv(pca_coords, file.path(out_dir, "Supplementary_Figure_S2_PCA_coordinates.csv"))
write_csv(qc_table, file.path(out_dir, "Supplementary_Table_S1_RNAseq_count_QC.csv"))

dose_colors <- c("0" = "#7F7F7F", "6.5" = "#6BAED6", "13" = "#3182BD", "19.5" = "#08519C")
gravity_shapes <- c("1g" = 16, "SMG" = 17)

plot_pca <- function(obj, label) {
  ggplot(obj$coords, aes(PC1, PC2, colour = dose_f, shape = gravity)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey75") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey75") +
    geom_point(size = 2.0, alpha = 0.9) +
    scale_colour_manual(values = dose_colors, name = "Dose (Gy)") +
    scale_shape_manual(values = gravity_shapes, name = "Gravity") +
    labs(
      x = sprintf("PC1 (%.1f%%)", obj$pct[1]),
      y = sprintf("PC2 (%.1f%%)", obj$pct[2]),
      title = label
    ) +
    theme_si() + theme(legend.position = "bottom")
}

s2 <- (plot_pca(pca3, "3 h") | plot_pca(pca10, "10 h")) +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_bundle(s2, "Supplementary_Figure_S2", 8.8)

## -------------------------------------------------------------------------
## Supplementary Figure S3 and Table S4: individual logistic fits
## -------------------------------------------------------------------------

od <- read_csv(file.path(data_root, "phenotype", "cleaned_OD600_data.csv"), show_col_types = FALSE) |>
  mutate(
    Dose = as.numeric(Dose),
    Dose_f = factor(Dose, levels = c(0, 6.5, 13, 19.5)),
    Condition = factor(Condition, levels = c("1g", "SMG"))
  )
fit_t <- c(0, 1, 2, 3, 4.5, 6, 7.5, 9, 10)
fit_cols <- c("0h", "1h", "2h", "3h", "4.5h", "6h", "7.5h", "9h", "10h")
logistic4 <- function(t, A, K, r, t0) A + (K - A) / (1 + exp(-r * (t - t0)))

fit_one <- function(i) {
  y <- as.numeric(od[i, fit_cols])
  tryCatch({
    fit <- nlsLM(
      y ~ logistic4(fit_t, A, K, r, t0),
      start = c(A = min(y), K = max(y), r = 0.8, t0 = 5),
      lower = c(A = 0, K = 0.05, r = 0.01, t0 = -5),
      upper = c(A = 0.2, K = 0.35, r = 3, t0 = 20),
      control = nls.lm.control(maxiter = 500)
    )
    cf <- coef(fit)
    pred <- logistic4(fit_t, cf["A"], cf["K"], cf["r"], cf["t0"])
    r2 <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)
    tibble(
      SampleID = od$SampleID[i], Condition = as.character(od$Condition[i]), Dose = od$Dose[i],
      A = unname(cf["A"]), K = unname(cf["K"]), r = unname(cf["r"]), t0 = unname(cf["t0"]),
      lag_time = unname(cf["t0"] - 2 / cf["r"]), R2 = r2, converged = TRUE
    )
  }, error = function(e) {
    tibble(
      SampleID = od$SampleID[i], Condition = as.character(od$Condition[i]), Dose = od$Dose[i],
      A = NA_real_, K = NA_real_, r = NA_real_, t0 = NA_real_, lag_time = NA_real_,
      R2 = NA_real_, converged = FALSE
    )
  })
}

fits <- bind_rows(lapply(seq_len(nrow(od)), fit_one)) |>
  mutate(
    Dose_f = factor(Dose, levels = c(0, 6.5, 13, 19.5)),
    Condition = factor(Condition, levels = c("1g", "SMG"))
  )
if (nrow(fits) != 24L || any(!fits$converged)) stop("Not all 24 logistic fits converged")
write_csv(
  fits |> mutate(across(c(A, K, r, t0, lag_time, R2), ~ round(.x, 4))),
  file.path(out_dir, "Supplementary_Table_S4_logistic_model_fit_quality.csv")
)

pred_t <- seq(0, 10, length.out = 201)
pred <- bind_rows(lapply(seq_len(nrow(fits)), function(i) {
  tibble(
    SampleID = fits$SampleID[i], Condition = fits$Condition[i], Dose_f = fits$Dose_f[i],
    Time = pred_t,
    OD_pred = logistic4(pred_t, fits$A[i], fits$K[i], fits$r[i], fits$t0[i])
  )
}))
raw_long <- od |>
  pivot_longer(all_of(fit_cols), names_to = "TimeLabel", values_to = "OD") |>
  mutate(Time = as.numeric(sub("h", "", TimeLabel)))

s3 <- ggplot() +
  geom_point(data = raw_long, aes(Time, OD), size = 0.55, alpha = 0.55, colour = "grey55") +
  geom_line(data = pred, aes(Time, OD_pred, group = SampleID, colour = Condition), linewidth = 0.45) +
  facet_grid(Condition ~ Dose_f, labeller = labeller(Dose_f = function(x) paste0(x, " Gy"))) +
  scale_colour_manual(values = c("1g" = "#1F78B4", "SMG" = "#D95F02")) +
  labs(x = "Time (h)", y = "OD₆₀₀") +
  theme_si() + theme(legend.position = "none")
save_bundle(s3, "Supplementary_Figure_S3", 11.8)

message("Supplementary Figures S1-S3 and Tables S1/S4 generated successfully.")
