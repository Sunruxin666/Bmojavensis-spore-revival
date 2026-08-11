# Figure 2: fixed microscopy montage plus reproducible quantitative panels
#
# The 32 raw microscopy fields used to assemble panel a are not available in
# this repository. Panel a is therefore stored from the author-approved,
# flattened manuscript figure and treated as a versioned visual input. Panels
# b and c are regenerated from biological-replicate germination measurements.

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(ggplot2)
  library(patchwork)
  library(png)
  library(readr)
  library(showtext)
  library(stringr)
  library(tidyr)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
input_file <- file.path(project_root, "data", "phenotype", "cleaned_germination_data.csv")
panel_a_file <- file.path(project_root, "data", "figure2", "Figure2A_microscopy_montage_flattened.png")
outdir <- file.path(project_root, "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

font_regular <- if (Sys.info()[["sysname"]] == "Darwin") {
  "/System/Library/Fonts/Supplemental/Arial.ttf"
} else {
  "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"
}
font_bold <- if (Sys.info()[["sysname"]] == "Darwin") {
  "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
} else {
  "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"
}
font_add("Arial", regular = font_regular, bold = font_bold)
showtext_auto()
showtext_opts(dpi = 600)

ci95 <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- mean(x)
  if (n <= 1) return(c(mean = m, lwr = NA, upr = NA, se = NA))
  se <- sd(x) / sqrt(n)
  tcrit <- qt(0.975, df = n - 1)
  c(mean = m, lwr = m - tcrit * se, upr = m + tcrit * se, se = se)
}

g <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    Dose = as.numeric(Dose),
    Dose_f = factor(Dose, levels = c(0, 6.5, 13, 19.5)),
    Condition = factor(Condition, levels = c("1g", "SMG"))
  )

traj <- g %>%
  transmute(
    SampleID, Condition, Dose, Dose_f,
    `0 h` = 0, `1 h` = Rate_1h, `2 h` = Rate_2h, `3 h` = Rate_3h
  ) %>%
  pivot_longer(
    cols = c(`0 h`, `1 h`, `2 h`, `3 h`),
    names_to = "TimeLabel", values_to = "Germ"
  ) %>%
  mutate(Time = as.numeric(str_remove(TimeLabel, " h")))

traj_sum <- traj %>%
  group_by(Condition, Dose, Dose_f, Time) %>%
  summarise(
    mean = ci95(Germ)[["mean"]],
    lwr = ci95(Germ)[["lwr"]],
    upr = ci95(Germ)[["upr"]],
    .groups = "drop"
  )

# Smithson-Verkuilen transformation keeps exact zeros inside (0, 1) for the
# beta mixed model while retaining SampleID as a repeated-measures intercept.
traj_beta <- traj %>%
  mutate(
    Nobs = n(),
    germ_beta = (Germ * (Nobs - 1) + 0.5) / Nobs
  )

fit <- glmmTMB(
  germ_beta ~ Time * Dose_f * Condition + (1 | SampleID),
  family = beta_family(link = "logit"), data = traj_beta
)

pred_grid <- expand.grid(
  Time = seq(0, 3, by = 0.05),
  Dose_f = levels(g$Dose_f),
  Condition = levels(g$Condition),
  SampleID = unique(g$SampleID)[1]
)
pred_link <- predict(
  fit, newdata = pred_grid, type = "link", se.fit = TRUE,
  re.form = NA, allow.new.levels = TRUE
)
pred_grid <- pred_grid %>%
  mutate(
    fit = plogis(pred_link$fit),
    lwr = plogis(pred_link$fit - 1.96 * pred_link$se.fit),
    upr = plogis(pred_link$fit + 1.96 * pred_link$se.fit)
  )

pred_diff <- pred_grid %>%
  select(Time, Dose_f, Condition, fit) %>%
  pivot_wider(names_from = Condition, values_from = fit) %>%
  mutate(diff = SMG - `1g`)

col_dose <- c("0" = "#0072B5", "6.5" = "#20854E", "13" = "#E18727", "19.5" = "#BC3C29")

theme_figure2 <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    axis.title = element_text(size = 8, colour = "black"),
    axis.text = element_text(size = 7, colour = "black"),
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks = element_line(linewidth = 0.4, colour = "black"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.width = unit(9, "pt"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    plot.tag = element_text(face = "bold", size = 9),
    plot.tag.position = c(0, 1),
    plot.margin = margin(4, 6, 4, 6)
  )

p_b <- ggplot() +
  geom_line(
    data = pred_grid,
    aes(Time, fit, colour = Dose_f, group = Dose_f),
    linewidth = 0.45
  ) +
  geom_point(
    data = traj_sum,
    aes(Time, mean, colour = Dose_f),
    size = 1.6
  ) +
  facet_wrap(~Condition, nrow = 1) +
  scale_colour_manual(values = col_dose) +
  scale_x_continuous(breaks = 0:3) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.8, 0.2)) +
  labs(tag = "b", x = "Sampling time [h]", y = "Relative germination rate") +
  theme_figure2

p_c <- ggplot(pred_diff, aes(Time, diff, colour = Dose_f, group = Dose_f)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "grey60") +
  geom_line(linewidth = 0.45) +
  scale_colour_manual(values = col_dose) +
  scale_x_continuous(breaks = 0:3) +
  labs(tag = "c", x = "Sampling time [h]", y = expression(Delta ~ "germination rate")) +
  theme_figure2

# The author-approved flattened panel-a input is 1259 x 535 px. Its 32 source
# microscopy fields are not available in this repository.
panel_a_raster <- png::readPNG(panel_a_file)
if (dim(panel_a_raster)[1] != 535L || dim(panel_a_raster)[2] != 1259L) {
  stop("Unexpected Figure 2a dimensions; expected 1259 x 535 px.")
}
p_a <- wrap_elements(full = grid::rasterGrob(panel_a_raster, interpolate = FALSE))

figure2 <- p_a / (p_b | p_c) +
  plot_layout(heights = c(0.82, 1), widths = c(1.05, 1)) &
  theme(plot.background = element_rect(fill = "white", colour = NA))

write_csv(traj_sum, file.path(outdir, "Figure2_observed_germination_summary.csv"))
write_csv(pred_grid, file.path(outdir, "Figure2_model_predictions.csv"))
write_csv(pred_diff, file.path(outdir, "Figure2_model_implied_SMG_minus_1g.csv"))

ggsave(
  file.path(outdir, "Figure2_submission.pdf"), figure2,
  width = 190, height = 180, units = "mm", device = cairo_pdf
)
ggsave(
  file.path(outdir, "Figure2_submission.png"), figure2,
  width = 190, height = 180, units = "mm", dpi = 600,
  device = ragg::agg_png, background = "white"
)

message("Figure 2 generated; panel a is a disclosed fixed raster input, panels b/c are data-derived.")
