project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

input_file <- file.path(
  project_root, "data", "phenotype", "cleaned_OD600_data.csv"
)
output_dir <- file.path(project_root, "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

od <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("SampleID", "Condition", "Dose", "0h")
missing_columns <- setdiff(required_columns, names(od))
if (length(missing_columns)) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}
if (anyNA(od[["0h"]]) || any(od[["0h"]] <= 0)) {
  stop("The 0 h OD600 values must be complete and positive.")
}

group_means <- aggregate(
  od[["0h"]],
  by = list(Condition = od$Condition, Dose_Gy = od$Dose),
  FUN = mean
)
names(group_means)[3] <- "mean_OD600_0h"
group_means <- group_means[order(group_means$Condition, group_means$Dose_Gy), ]

summary_metrics <- data.frame(
  metric = c(
    "overall_mean_OD600_0h",
    "CV_percent_among_eight_treatment_group_means",
    "CV_percent_among_24_individual_values"
  ),
  value = c(
    mean(od[["0h"]]),
    100 * sd(group_means$mean_OD600_0h) / mean(group_means$mean_OD600_0h),
    100 * sd(od[["0h"]]) / mean(od[["0h"]])
  ),
  n = c(nrow(od), nrow(group_means), nrow(od))
)

write.csv(
  group_means,
  file.path(output_dir, "inoculum_OD600_treatment_group_means.csv"),
  row.names = FALSE
)
write.csv(
  summary_metrics,
  file.path(output_dir, "inoculum_OD600_uniformity_qc.csv"),
  row.names = FALSE
)

group_cv <- summary_metrics$value[
  summary_metrics$metric == "CV_percent_among_eight_treatment_group_means"
]
message(
  sprintf(
    "Inoculum OD600 QC: overall mean = %.6f; CV among eight treatment-group means = %.3f%%.",
    mean(od[["0h"]]), group_cv
  )
)

