#!/usr/bin/env Rscript

# Assemble the two data-derived Figure 6 panels after their individual scripts
# have run. The former interpretive schematic panel has been removed.

suppressPackageStartupMessages({
  library(png)
  library(ragg)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
out_dir <- file.path(project_root, "results")

panel_files <- c(
  a = file.path(out_dir, "Figure6a_submission_v16.png"),
  b = file.path(out_dir, "Figure6b_submission_v16.png")
)

missing <- panel_files[!file.exists(panel_files)]
if (length(missing)) {
  stop("Missing Figure 6 panel exports:\n- ", paste(missing, collapse = "\n- "))
}

panels <- lapply(panel_files, png::readPNG)
target_width <- dim(panels$a)[2]
if (any(vapply(panels, function(x) dim(x)[2], integer(1)) != target_width)) {
  stop("Figure 6 panels must share the same pixel width before assembly.")
}

draw_integrated_figure6 <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 2, ncol = 1,
    heights = grid::unit(c(14, 16), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))

  for (i in seq_along(panels)) {
    grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
    grid::grid.raster(panels[[i]], width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
                      interpolate = FALSE)
    grid::grid.text(
      names(panels)[i], x = grid::unit(0.012, "npc"), y = grid::unit(0.985, "npc"),
      just = c("left", "top"), gp = grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9)
    )
    grid::popViewport()
  }
  grid::popViewport()
}

ragg::agg_png(
  file.path(out_dir, "Figure6_integrated_submission.png"),
  width = 18.3, height = 30, units = "cm", res = 600, background = "white"
)
draw_integrated_figure6()
grDevices::dev.off()

ragg::agg_tiff(
  file.path(out_dir, "Figure6_integrated_submission.tiff"),
  width = 18.3, height = 30, units = "cm", res = 600,
  compression = "lzw", background = "white"
)
draw_integrated_figure6()
grDevices::dev.off()

grDevices::cairo_pdf(
  file.path(out_dir, "Figure6_integrated_submission.pdf"),
  width = 18.3 / 2.54, height = 30 / 2.54, onefile = TRUE
)
draw_integrated_figure6()
grDevices::dev.off()

ragg::agg_jpeg(
  file.path(out_dir, "Figure6_integrated_submission.jpg"),
  width = 18.3, height = 30, units = "cm", res = 600,
  quality = 97, background = "white"
)
draw_integrated_figure6()
grDevices::dev.off()

message("Integrated Figure 6 generated from data-derived panels a-b.")
