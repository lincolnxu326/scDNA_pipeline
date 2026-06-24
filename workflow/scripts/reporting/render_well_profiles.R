#!/usr/bin/env Rscript

# Render one AneuFinder copy-number profile PNG per well for the QC review report.
#
# This is a read-only consumer of the per-well model objects already written by
# run_aneufinder.R. It does NOT re-run AneuFinder or touch its outputs; it only
# loads each MODELS/method-{method}/{well}.RData object and re-renders the same
# profile plot used in profiles_{method}.pdf as a single PNG per well, plus a
# manifest.json the HTML generator consumes. Keep it decoupled from
# run_aneufinder.R (which is fragile / mid-debug).

if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

library(AneuFinder)
library(optparse)
library(logger)

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "AneuFinder MODELS directory (contains method-{method}/)"),
  make_option(c("-o", "--outdir"), type = "character",
              help = "Output directory for per-well PNGs and manifest.json"),
  make_option(c("-m", "--method"), type = "character", default = "edivisive",
              help = "Segmentation method subdirectory to render [default %default]"),
  make_option(c("--heatmap"), type = "character", default = NULL,
              help = "Optional path; if set, also render a genome-wide CN heatmap PNG"),
  make_option(c("--format"), type = "character", default = "png",
              help = "Per-well plot format: 'png' (compact) or 'svg' (crisp/vector) [default %default]"),
  make_option(c("--width"), type = "integer", default = 1400,
              help = "PNG width in pixels [default %default]"),
  make_option(c("--height"), type = "integer", default = 430,
              help = "PNG height in pixels [default %default]"),
  make_option(c("--res"), type = "integer", default = 110,
              help = "PNG resolution in ppi [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Minimal JSON string escaper so we avoid depending on jsonlite in the renv path.
json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  x
}

# Write a simple placeholder PNG carrying a message, so a failed heatmap never
# breaks the downstream HTML viewer.
write_placeholder_png <- function(path, message, width = 1200, height = 400, res = 110) {
  tryCatch({
    grDevices::png(filename = path, width = width, height = height, res = res)
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::plot.new()
    graphics::text(0.5, 0.5, message, cex = 1.2)
    grDevices::dev.off()
  }, error = function(e) {
    try(grDevices::dev.off(), silent = TRUE)
  })
}

# Render a genome-wide copy-number heatmap across all rendered models to a PNG.
# heatmapGenomewide() accepts a vector of model filenames and returns a ggplot.
render_genome_heatmap <- function(rdata_files, path) {
  if (length(rdata_files) == 0) {
    write_placeholder_png(path, "No models available for genome heatmap")
    return(invisible(FALSE))
  }
  hp <- tryCatch(
    AneuFinder::heatmapGenomewide(rdata_files),
    error = function(err) {
      log_warn(paste0("Could not build genome heatmap: ", err$message))
      NULL
    }
  )
  if (is.null(hp)) {
    write_placeholder_png(path, "Genome heatmap unavailable")
    return(invisible(FALSE))
  }
  ok <- tryCatch({
    # very wide so per-chromosome columns + labels are legible across the genome
    grDevices::png(filename = path, width = 3000,
                   height = max(420, 42 * length(rdata_files)), res = 110)
    print(hp)
    grDevices::dev.off()
    TRUE
  }, error = function(err) {
    try(grDevices::dev.off(), silent = TRUE)
    log_warn(paste0("Could not write genome heatmap PNG: ", err$message))
    FALSE
  })
  if (!ok) write_placeholder_png(path, "Genome heatmap unavailable")
  invisible(ok)
}

render_well_profiles <- function() {
  if (is.null(opt$input))  stop("--input (MODELS directory) is required")
  if (is.null(opt$outdir)) stop("--outdir is required")

  method_dir <- file.path(opt$input, paste0("method-", opt$method))
  if (!dir.exists(method_dir)) {
    stop(paste0("Method directory not found: ", method_dir))
  }

  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

  rdata_files <- list.files(method_dir, pattern = "\\.RData$", full.names = TRUE)
  rdata_files <- sort(rdata_files)
  if (length(rdata_files) == 0) {
    log_warn(paste0("No .RData models found in ", method_dir))
  }
  log_info(paste0("Rendering ", length(rdata_files), " profile plots from ", method_dir))

  # Render an AneuFinder plot of the given type, wide aspect for a readable genome
  # x-axis. format="svg" → vector (crisp/responsive; large files — good for the small
  # cn_review). format="png" → raster (compact; good for the 96-well first-pass review).
  # SVG always falls back to PNG if no SVG device is available. Returns basename or NA.
  PNG_DPI <- 150L   # wide, crisp raster
  render_vec <- function(model, type, base_noext, w_in, h_in, breakpoints = FALSE,
                         format = "png") {
    p <- tryCatch(
      if (breakpoints) plot(model, type = type, plot.breakpoints = TRUE)
      else plot(model, type = type),
      error = function(err) {
        log_warn(paste0("plot(", type, ") failed for ", basename(base_noext),
                        ": ", err$message)); NULL
      })
    if (is.null(p)) return(NA_character_)
    render_png <- function() {
      png_path <- paste0(base_noext, ".png")
      ok <- tryCatch({
        grDevices::png(filename = png_path, width = round(w_in * PNG_DPI),
                       height = round(h_in * PNG_DPI), res = PNG_DPI)
        print(p); grDevices::dev.off(); TRUE
      }, error = function(err) {
        try(grDevices::dev.off(), silent = TRUE)
        log_warn(paste0("Could not render ", type, " for ", basename(base_noext),
                        ": ", err$message)); FALSE
      })
      if (ok) basename(png_path) else NA_character_
    }
    if (format == "svg") {
      svg_path <- paste0(base_noext, ".svg")
      ok_svg <- tryCatch({
        if (requireNamespace("svglite", quietly = TRUE)) {
          svglite::svglite(filename = svg_path, width = w_in, height = h_in)
        } else {
          grDevices::svg(filename = svg_path, width = w_in, height = h_in)
        }
        print(p); grDevices::dev.off(); TRUE
      }, error = function(err) { try(grDevices::dev.off(), silent = TRUE); FALSE })
      if (ok_svg && file.exists(svg_path) && file.info(svg_path)$size > 0)
        return(basename(svg_path))
      if (file.exists(svg_path)) unlink(svg_path)   # drop empty/partial svg
      log_warn(paste0("SVG unavailable for ", basename(base_noext), "; using PNG"))
    }
    render_png()
  }

  entries <- list()
  n_ok <- 0
  for (ifile in rdata_files) {
    well <- tools::file_path_sans_ext(basename(ifile))
    profile_file <- NA_character_
    hist_file <- NA_character_
    tryCatch({
      obj_name <- load(ifile)
      model <- get(obj_name[1])
      # profile (same call as run_aneufinder.R profiles PDF) — wide for a readable x-axis
      # very wide aspect so the genome x-axis (chromosome labels) is not crammed
      profile_file <- render_vec(model, "profile",
                                 file.path(opt$outdir, paste0(well, "_profile")),
                                 24, 3.4, breakpoints = TRUE, format = opt$format)
      # bin read-count histogram with fitted somy/state densities
      hist_file <- render_vec(model, "histogram",
                              file.path(opt$outdir, paste0(well, "_histogram")), 8, 3.4,
                              format = opt$format)
    }, error = function(err) {
      try(grDevices::dev.off(), silent = TRUE)
      log_warn(paste0("Could not load model for ", well, ": ", err$message))
    })
    if (!is.na(profile_file)) n_ok <- n_ok + 1
    j <- function(x) if (is.na(x)) "null" else paste0("\"", json_escape(x), "\"")
    entries[[length(entries) + 1]] <- paste0(
      "{\"well\": \"", json_escape(well),
      "\", \"profile\": ", j(profile_file),
      ", \"histogram\": ", j(hist_file),
      ", \"ok_profile\": ", if (!is.na(profile_file)) "true" else "false",
      ", \"ok_histogram\": ", if (!is.na(hist_file)) "true" else "false", "}"
    )
  }

  manifest_path <- file.path(opt$outdir, "manifest.json")
  manifest <- paste0(
    "{\n",
    "  \"method\": \"", json_escape(opt$method), "\",\n",
    "  \"wells\": [\n    ",
    paste(entries, collapse = ",\n    "),
    "\n  ]\n}\n"
  )
  writeLines(manifest, manifest_path)
  log_info(paste0("Wrote manifest: ", manifest_path))
  log_info(paste0("Rendered ", n_ok, "/", length(rdata_files),
                  " per-well profile plots to ", opt$outdir))

  # Optional genome-wide CN heatmap (used by the second-pass / cn_review report).
  if (!is.null(opt$heatmap)) {
    dir.create(dirname(opt$heatmap), recursive = TRUE, showWarnings = FALSE)
    render_genome_heatmap(rdata_files, opt$heatmap)
    log_info(paste0("Genome heatmap written to: ", opt$heatmap))
  }
}

if (!interactive()) {
  tryCatch(render_well_profiles(), error = function(e) {
    log_error(e$message)
    quit(status = 1)
  })
}
