project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(project_root, "R", "preflight.R"), local = new.env(parent = globalenv()))

scripts <- c(
  "figure2_submission.R",
  "figure3_submission.R",
  "figure4_submission.R",
  "figure5_submission.R",
  "figure6a_submission.R",
  "figure6b_submission.R",
  "figure6_integrate.R",
  "figure7_submission.R",
  "supplementary_figures_submission.R",
  "figure34_statistical_audit.R",
  "figure7_replicate_auc_audit.R"
)

for (script in scripts) {
  message("\n=== Running ", script, " ===")
  source(
    file.path(project_root, "R", script),
    local = new.env(parent = globalenv()),
    chdir = FALSE
  )
}

writeLines(
  capture.output(sessionInfo()),
  file.path(project_root, "results", "sessionInfo.txt")
)
message("\nAll figure and audit scripts completed.")
