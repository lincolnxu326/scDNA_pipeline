#!/usr/bin/env Rscript

## Activate renv if present
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

suppressPackageStartupMessages({
  library(AneuFinder)
  library(optparse)
  library(logger)
  library(GenomicRanges)
})

## Fixed path to the variable-width reference BED
VARIABLE_WIDTH_REFERENCE <- "/nemo/project/proj-tracerX/working/VCAM1_GnT/DAQI/source/Aneufinder/mappability/normal_ref_sorted.bam.bed"

## CLI
option_list <- list(
  make_option(c("-i", "--input"), type = "character"),
  make_option(c("-o", "--output"), type = "character"),
  make_option(c("-b", "--blacklist"), type = "character", default = NULL),
  make_option(c("-m", "--method"), type = "character", default = "edivisive"),
  make_option(c("-s", "--binsize"), type = "integer", default = 1e6),
  make_option(c("-c", "--chromosomes"), type = "character", default = paste0("chr", c(1:22))),
  make_option(c("-n", "--numcpu"), type = "integer", default = 4),
  make_option(c("-p", "--cluster-plots"), type = "logical", default = TRUE),
  make_option(c("-r", "--reuse-existing"), type = "logical", default = FALSE),
  # GC template is now REQUIRED
  make_option(c("-g", "--gc-rds"), type = "character", default = NULL,
              help = "Plate-specific GC template RDS (generated in relaxed env). REQUIRED.")
)
opt <- parse_args(OptionParser(option_list = option_list))

run_aneufinder <- function() {
  
  # hard fail if GC not provided
  if (is.null(opt$`gc-rds`)) {
    stop("You must provide --gc-rds pointing to the plate-specific GC template .rds")
  }
  if (!file.exists(opt$`gc-rds`)) {
    stop(paste0("GC template file not found: ", opt$`gc-rds`))
  }
  
  log_info("Starting AneuFinder analysis")
  log_info(paste0("Input: ", opt$input))
  log_info(paste0("Output: ", opt$output))
  log_info(paste0("GC template: ", opt$`gc-rds`))
  
  if (!dir.exists(opt$input)) stop("Input directory not found: ", opt$input)
  
  bam_files <- list.files(opt$input, pattern = "\\.bam$", full.names = TRUE)
  if (length(bam_files) == 0) stop("No BAM files found in: ", opt$input)
  
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  
  chromosomes <- unlist(strsplit(opt$chromosomes, ","))
  methods     <- unlist(strsplit(opt$method, ","))
  
  ## 1) Precompute variable width bins once from the big BED
  log_info("Precomputing variable width bins once from reference BED ...")
  ref_reads_ <- bed2GRanges(
    bedfile    = VARIABLE_WIDTH_REFERENCE,
    assembly   = "hg38",
    chromosomes = chromosomes,
    blacklist  = opt$blacklist
  )
  
  pre_bins_ <- variableWidthBins(
    reads       = ref_reads_,
    binsizes    = opt$binsize,
    chromosomes = chromosomes
  )
  
  ## 2) Load plate-specific GC template (must exist)
  bins_gc_ <- readRDS(opt$`gc-rds`)
  if (inherits(bins_gc_, "GRanges")) {
    bins_gc_ <- GenomicRanges::GRangesList("0" = bins_gc_)
  }
  
  # make model dirs up front
  for (m in methods) {
    dir.create(file.path(opt$output, "MODELS", paste0("method-", m)),
               recursive = TRUE, showWarnings = FALSE)
  }
  
  ## 3) Process BAMs
  for (bf in bam_files) {
    id_ <- tools::file_path_sans_ext(basename(bf))
    log_info(paste0("Processing ", id_))
    
    # bin using precomputed bins
    binned_ <- binReads(
      file           = bf,
      bins           = pre_bins_,
      blacklist      = opt$blacklist,
      assembly       = "hg38",
      save.as.RData  = FALSE,
      use.bamsignals = FALSE
    )
    
    # normalise to GRangesList
    if (inherits(binned_, "GRanges") || inherits(binned_, "GRangesList")) {
      binned_obj_ <- binned_
    } else {
      binned_obj_ <- binned_[[1]]
    }
    if (inherits(binned_obj_, "GRanges")) {
      binned_obj_ <- GenomicRanges::GRangesList("0" = binned_obj_)
    }
    
    # GC correction with the plate-specific template
    corrected_list_ <- correctGC(
      binned.data.list = list(binned_obj_),
      GC.BSgenome      = NULL,
      same.binsize     = TRUE,
      method           = "loess",
      bins             = bins_gc_
    )
    binned_gc_obj_ <- corrected_list_[[1]]
    
    # CNV calling
    for (m in methods) {
      method_dir_ <- file.path(opt$output, "MODELS", paste0("method-", m))
      dir.create(method_dir_, recursive = TRUE, showWarnings = FALSE)
      
      model <- switch(
        m,
        "edivisive" = findCNVs.strandseq(binned_gc_obj_, ID = id_, method = "edivisive"),
        "HMM"       = findCNVs.strandseq(binned_gc_obj_, ID = id_, method = "HMM"),
        "dnacopy"   = findCNVs.strandseq(binned_gc_obj_, ID = id_, method = "dnacopy"),
        stop(paste("Unknown method:", m))
      )
      save(model, file = file.path(method_dir_, paste0(id_, ".RData")))
    }
  }
  
  ## 4) Plotting (your block)
  if (isTRUE(opt$`cluster-plots`)) {
    for (m in methods) {
      method_dir <- file.path(opt$output, "MODELS", paste0("method-", m))
      if (!dir.exists(method_dir)) next
      
      rdata_files <- list.files(method_dir, pattern = "\\.RData$", full.names = TRUE)
      if (length(rdata_files) == 0) next
      
      clus <- clusterByQuality(rdata_files)
      out_pdf <- file.path(opt$output, paste0("Genome_heatmap_cluster_1Mb_bins_", m, ".pdf"))
      heatmapGenomewideClusters(clus, file = out_pdf)
      log_info(paste0("Cluster heatmap written to: ", out_pdf))
    }
  }
  
  ## 5) Summary
  results_dir <- file.path(opt$output, "MODELS")
  if (dir.exists(results_dir)) {
    model_count  <- length(list.files(results_dir, pattern = "\\.RData$", recursive = TRUE))
    summary_file <- file.path(opt$output, "analysis_summary.txt")
    writeLines(c(
      "ANEUFINDER ANALYSIS SUMMARY",
      "===========================",
      paste0("Date: ", Sys.Date()),
      paste0("Input directory: ", opt$input),
      paste0("Total BAM files: ", length(bam_files)),
      paste0("Successful analyses: ", model_count),
      paste0("Methods used: ", paste(methods, collapse = ", ")),
      paste0("Bin size: ", format(opt$binsize, big.mark = ",")),
      paste0("Chromosomes: ", paste(chromosomes, collapse = ", ")),
      if (!is.null(opt$blacklist)) paste0("Blacklist: ", opt$blacklist) else NULL,
      paste0("Variable width reference: ", VARIABLE_WIDTH_REFERENCE),
      paste0("GC template: ", opt$`gc-rds`)
    ), summary_file)
    log_info(paste0("Summary written to: ", summary_file))
  }
}

if (!interactive()) {
  tryCatch(run_aneufinder(), error = function(e) {
    log_error(e$message)
    quit(status = 1)
  })
}
