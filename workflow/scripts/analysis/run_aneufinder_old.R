#!/usr/bin/env Rscript
#' Run AneuFinder for single-cell copy number analysis
#'
#' This script runs AneuFinder on aligned BAM files to detect
#' copy number variations in single cells.
#'
#' NOTE: GC correction is currently non-functional due to BSgenome.Hsapiens.UCSC.hg38
#' not being available in the current environment.

# Activate renv for reproducible R environment
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

# Set download method for restricted networks
options(download.file.method = "wget")

suppressPackageStartupMessages({
  # Load AneuFinder - should be available through renv
  if (!requireNamespace("AneuFinder", quietly = TRUE)) {
    stop("AneuFinder is required but not found in renv library. Please check renv installation.")
  }
  library(AneuFinder)

  # Note: BSgenome.Hsapiens.UCSC.hg38 is NOT available
  # GC correction will be disabled

  library(optparse)
  library(logger)
})

# Variable bin width, suitable for plate 7. 
VARIABLE_WIDTH_REFERENCE <- "/nemo/project/proj-tracerX/working/VCAM1_GnT/DAQI/source/Aneufinder/mappability/normal_ref_sorted.bam.bed"

# Parse command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Input directory containing BAM files"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output directory for AneuFinder results"),
  make_option(c("-b", "--blacklist"), type = "character", default = NULL,
              help = "Path to blacklist BED file"),
  make_option(c("-m", "--method"), type = "character", default = "edivisive",
              help = "Method(s) for segmentation (comma-separated)"),
  make_option(c("-s", "--binsize"), type = "integer", default = 1e6,
              help = "Bin size for analysis"),
  make_option(c("-c", "--chromosomes"), type = "character",
              default = paste0("chr", c(1:22)),
              help = "Chromosomes to analyze (comma-separated)"),
  make_option(c("-n", "--numcpu"), type = "integer", default = 4,
              help = "Number of CPUs to use"),
  make_option(c("-p", "--cluster-plots"), type = "logical", default = TRUE,
              help = "Generate cluster plots"),
  make_option(c("-r", "--reuse-existing"), type = "logical", default = FALSE,
              help = "Reuse existing results if available")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Main function
run_aneufinder <- function() {
  log_info("Starting AneuFinder analysis")
  log_info(paste0("Input directory: ", opt$input))
  log_info(paste0("Output directory: ", opt$output))

  # Check input directory
  if (!dir.exists(opt$input)) {
    log_error(paste0("Input directory does not exist: ", opt$input))
    stop("Input directory not found")
  }

  # List BAM files
  bam_files <- list.files(opt$input, pattern = "\\.bam$", full.names = TRUE)
  log_info(paste0("Found ", length(bam_files), " BAM files"))

  if (length(bam_files) == 0) {
    log_error("No BAM files found in input directory")
    stop("No BAM files found")
  }

  # Create output directory
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)

  # Parse chromosomes
  chromosomes <- unlist(strsplit(opt$chromosomes, ","))

  # Parse methods
  methods <- unlist(strsplit(opt$method, ","))

  # Prepare AneuFinder parameters
  params <- list(
    inputfolder = opt$input,
    outputfolder = opt$output,
    numCPU = opt$numcpu,
    method = methods,
    binsizes = opt$binsize,
    chromosomes = chromosomes,
    reuse.existing.files = opt$`reuse-existing`
  )

  # Add blacklist if provided
  if (!is.null(opt$blacklist) && file.exists(opt$blacklist)) {
    log_info(paste0("Using blacklist: ", opt$blacklist))
    params$blacklist <- opt$blacklist
  }

  # Add GC correction parameters if BSgenome was available (currently disabled)
  log_warn("GC correction is disabled - BSgenome.Hsapiens.UCSC.hg38 not available")

  # Add variable width reference (your fixed path)
  if (file.exists(VARIABLE_WIDTH_REFERENCE)) {
    log_info(paste0("Using variable.width.reference: ", VARIABLE_WIDTH_REFERENCE))
    params$variable.width.reference <- VARIABLE_WIDTH_REFERENCE
  } else {
    log_warn(paste0("variable.width.reference file not found: ", VARIABLE_WIDTH_REFERENCE,
                    " -- running without it"))
  }

  # Run AneuFinder
  log_info("Running AneuFinder with parameters:")
  log_info(paste0("  Methods: ", paste(methods, collapse = ", ")))
  log_info(paste0("  Bin size: ", opt$binsize))
  log_info(paste0("  Chromosomes: ", paste(chromosomes, collapse = ", ")))
  log_info(paste0("  CPUs: ", opt$numcpu))

  tryCatch({
    # Run AneuFinder
    do.call(Aneufinder, params)

    log_info("AneuFinder analysis completed successfully")

    # Generate additional plots if requested
    if (opt$`cluster-plots`) {
      log_info("Generating cluster plots")

      # Load results
      results_dir <- file.path(opt$output, "MODELS", paste0("method-", methods[1]))

      if (dir.exists(results_dir)) {
        # Get list of RData files
        model_files <- list.files(results_dir, pattern = "\\.RData$",
                                  full.names = TRUE, recursive = TRUE)

        if (length(model_files) > 0) {
          # Load models for clustering
          models <- list()
          for (i in seq_along(model_files)) {
            load(model_files[i])
            # The loaded object should be named 'model'
            if (exists("model")) {
              models[[basename(model_files[i])]] <- model
            }
          }

          if (length(models) > 1) {
            # Generate heatmap
            pdf(file.path(opt$output, "cluster_heatmap.pdf"),
                width = 12, height = 8)

            tryCatch({
              heatmapGenomewide(models, cluster = TRUE)
            }, error = function(e) {
              log_warn(paste0("Could not generate heatmap: ", e$message))
            })

            dev.off()

            log_info("Cluster plots generated")
          } else {
            log_warn("Not enough models for clustering")
          }
        }
      }
    }

  }, error = function(e) {
    log_error(paste0("AneuFinder failed: ", e$message))
    stop(e)
  })

  # Summary statistics
  log_info("Generating summary statistics")

  # Count successful analyses
  results_dir <- file.path(opt$output, "MODELS")
  if (dir.exists(results_dir)) {
    model_count <- length(list.files(results_dir, pattern = "\\.RData$",
                                     recursive = TRUE))
    log_info(paste0("Successfully analyzed ", model_count, " cells"))

    # Write summary
    summary_file <- file.path(opt$output, "analysis_summary.txt")
    summary_lines <- c(
      "ANEUFINDER ANALYSIS SUMMARY",
      "===========================",
      paste0("Date: ", Sys.Date()),
      paste0("Input directory: ", opt$input),
      paste0("Total BAM files: ", length(bam_files)),
      paste0("Successful analyses: ", model_count),
      paste0("Methods used: ", paste(methods, collapse = ", ")),
      paste0("Bin size: ", format(opt$binsize, big.mark = ",")),
      paste0("Chromosomes: ", paste(chromosomes, collapse = ", ")),
      "",
      "NOTES:",
      "- GC correction is currently disabled (BSgenome.Hsapiens.UCSC.hg38 not available)",
      "- This may affect copy number calling accuracy",
      ""
    )

    if (!is.null(opt$blacklist)) {
      summary_lines <- c(summary_lines,
                         paste0("Blacklist used: ", opt$blacklist))
    }

    if (file.exists(VARIABLE_WIDTH_REFERENCE)) {
      summary_lines <- c(summary_lines,
                         paste0("Variable width reference: ", VARIABLE_WIDTH_REFERENCE))
    }

    writeLines(summary_lines, summary_file)
    log_info(paste0("Summary written to: ", summary_file))
  }
}

# Run main function
if (!interactive()) {
  tryCatch({
    run_aneufinder()
  }, error = function(e) {
    log_error(paste0("Fatal error: ", e$message))
    quit(status = 1)
  })
}
