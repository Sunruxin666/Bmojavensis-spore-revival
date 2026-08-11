packages <- c(
  "car", "cowplot", "dplyr", "emmeans", "forcats", "ggalluvial",
  "ggplot2", "glmmTMB", "minpack.lm", "openxlsx", "patchwork", "png", "ragg", "readr",
  "readxl", "reshape2", "scales", "showtext", "stringr", "svglite", "tidyr"
)

missing <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required packages are already installed.")
}
