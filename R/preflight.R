project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)

required_packages <- c(
  "car", "cowplot", "dplyr", "emmeans", "forcats", "ggalluvial",
  "ggplot2", "glmmTMB", "minpack.lm", "openxlsx", "patchwork", "png", "ragg", "readr",
  "readxl", "reshape2", "scales", "showtext", "stringr", "svglite", "tidyr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages)) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nRun: Rscript requirements.R"
  )
}

required_files <- c(
  "data/figure2/Figure2A_microscopy_montage_flattened.png",
  "data/phenotype/cleaned_germination_data.csv",
  "data/phenotype/cleaned_OD600_data.csv",
  "data/phenotype/Fig5A_permanova_3h.csv",
  "data/phenotype/Fig5A_permanova_10h.csv",
  "data/phenotype/growth_kinetics_group_summary.csv",
  "data/annotations/ko_definition.tsv",
  "data/annotations/ko_to_pathway.tsv",
  "data/annotations/pathway_names.tsv",
  file.path("data/secondary_spores", c(
    "control_1.xlsx", "control_2.xlsx", "h2o2.xlsx", "gentamicin.xlsx",
    "DDAC.xlsx", "soda_1.xlsx"
  ))
)

missing_files <- required_files[!file.exists(file.path(project_root, required_files))]
if (length(missing_files)) {
  stop("Missing required inputs:\n- ", paste(missing_files, collapse = "\n- "))
}

deg_files <- list.files(
  file.path(project_root, "data", "transcriptomics"),
  pattern = "DEG_full_.*_vs_Ctrl[.]csv$", recursive = TRUE, full.names = TRUE
)
kegg3_files <- list.files(
  file.path(project_root, "data", "transcriptomics"),
  pattern = "all_genes[.]kegg3$", recursive = TRUE, full.names = TRUE
)

if (length(deg_files) != 14L) {
  stop("Expected 14 DEG contrast files, found ", length(deg_files))
}
if (length(kegg3_files) != 6L) {
  stop("Expected 6 gene-to-KO mapping files, found ", length(kegg3_files))
}

dir.create(file.path(project_root, "results"), recursive = TRUE, showWarnings = FALSE)
message("Preflight passed: packages, core inputs, 14 DEG files and 6 KO maps found.")
