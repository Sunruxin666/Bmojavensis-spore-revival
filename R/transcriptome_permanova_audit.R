#!/usr/bin/env Rscript

# Reproduce the stage-stratified transcriptome PERMANOVA and test homogeneity
# of multivariate dispersion from the count matrices distributed with this
# repository. The analysis follows the original low-count filter and Hellinger
# transformation used for Figure 5a.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(vegan)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
count_root <- file.path(project_root, "data", "transcriptomics", "count_matrices")
out_dir <- file.path(project_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

min_rowsum <- 20
permutations <- 999
dose_levels <- c(0, 6.5, 13, 19.5)

merge_count_files <- function(timepoint) {
  paths <- file.path(
    count_root, timepoint,
    c("Ctrl_all_CountMatrix.tsv", "Ctrl_vs_SMG_CountMatrix.tsv", "Ctrl_vs_SMG_SR_CountMatrix.tsv")
  )
  if (!all(file.exists(paths))) stop("Missing count matrix for ", timepoint)
  mats <- lapply(paths, read_tsv, show_col_types = FALSE)
  genes <- Reduce(intersect, lapply(mats, function(x) x$Geneid))
  mats <- lapply(mats, function(x) as.data.frame(x[match(genes, x$Geneid), ]))

  counts <- data.frame(Geneid = genes, check.names = FALSE)
  sample_names <- unique(unlist(lapply(mats, function(x) setdiff(names(x), "Geneid"))))
  for (sample_name in sample_names) {
    candidates <- lapply(mats, function(x) if (sample_name %in% names(x)) x[[sample_name]] else NULL)
    candidates <- candidates[!vapply(candidates, is.null, logical(1))]
    if (length(candidates) > 1L &&
        !all(vapply(candidates[-1], identical, logical(1), candidates[[1]]))) {
      stop("Conflicting duplicated count column: ", sample_name, " at ", timepoint)
    }
    counts[[sample_name]] <- as.numeric(candidates[[1]])
  }

  metadata <- tibble(sample = sample_names) %>%
    mutate(
      condition = factor(if_else(str_detect(sample, "^SMG"), "SMG", "1g"),
                         levels = c("1g", "SMG")),
      dose = case_when(
        str_detect(sample, "19_5") ~ 19.5,
        str_detect(sample, "6_5") ~ 6.5,
        str_detect(sample, "13") ~ 13,
        TRUE ~ 0
      ),
      dose_f = factor(dose, levels = dose_levels)
    )

  matrix_counts <- as.matrix(counts[, metadata$sample, drop = FALSE])
  storage.mode(matrix_counts) <- "numeric"
  matrix_counts <- matrix_counts[rowSums(matrix_counts) >= min_rowsum, , drop = FALSE]
  transformed <- decostand(t(matrix_counts), method = "hellinger")
  stopifnot(identical(rownames(transformed), metadata$sample))
  list(expression = transformed, metadata = metadata)
}

run_permanova <- function(object) {
  # The RDA term test is retained here because it preceded PERMANOVA in the
  # original analysis and therefore preserves the published permutation stream.
  rda_fit <- rda(object$expression ~ condition * dose_f, data = object$metadata)
  invisible(anova.cca(rda_fit, by = "terms", permutations = permutations))
  distance <- vegdist(object$expression, method = "bray")
  fit <- adonis2(
    distance ~ condition * dose_f,
    data = object$metadata,
    by = "terms",
    permutations = permutations
  )
  list(distance = distance, fit = fit)
}

as_result_table <- function(x, timepoint) {
  as.data.frame(x) %>%
    rownames_to_column("term") %>%
    mutate(timepoint = timepoint, .before = 1)
}

objects <- list(`3 h` = merge_count_files("3h"), `10 h` = merge_count_files("10h"))

set.seed(123)
permanova_results <- lapply(objects, run_permanova)
permanova_table <- bind_rows(Map(
  function(result, label) as_result_table(result$fit, label),
  permanova_results, names(permanova_results)
))

# Use a separate declared seed so the diagnostic test does not alter the
# published PERMANOVA permutation stream.
set.seed(20260905)
permdisp_table <- bind_rows(Map(function(object, label) {
  distance <- vegdist(object$expression, method = "bray")
  groups <- interaction(object$metadata$condition, object$metadata$dose_f, drop = TRUE)
  dispersion_fit <- betadisper(distance, groups)
  permutation_test <- permutest(dispersion_fit, permutations = permutations)
  as.data.frame(permutation_test$tab) %>%
    rownames_to_column("term") %>%
    filter(term == "Groups") %>%
    transmute(
      timepoint = label,
      groups = nlevels(groups),
      Df,
      SumOfSqs = `Sum Sq`,
      MeanSqs = `Mean Sq`,
      F,
      permutations = N.Perm,
      p_value = `Pr(>F)`
    )
}, objects, names(objects)))

write_csv(permanova_table, file.path(out_dir, "Figure5_PERMANOVA_reproduced.csv"))
write_csv(permdisp_table, file.path(out_dir, "Figure5_PERMDISP_audit.csv"))

message("Transcriptome PERMANOVA and PERMDISP audit completed.")
