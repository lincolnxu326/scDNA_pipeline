#!/usr/bin/env Rscript
#' Filter AneuFinder results based on chr3p neutrality
#' 
#' This script filters cells based on chr3p copy number status
#' to remove potentially problematic cells from downstream analysis

# Activate renv for reproducible R environment
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

suppressPackageStartupMessages({
  library(AneuFinder)
  library(GenomicRanges)
  library(tidyverse)
  library(optparse)
  library(logger)
})

# Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "AneuFinder output directory"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output directory for filtered results"),
  make_option(c("-a", "--armlevel"), type = "character",
              default = "hg38.armlevel.cytoBand.txt",
              help = "Arm-level cytoband file"),
  make_option(c("-p", "--prop-threshold"), type = "numeric", default = 0.2,
              help = "Proportion threshold for segment consideration"),
  make_option(c("-c", "--chromosome"), type = "character", default = "chr3p",
              help = "Chromosome arm to check for neutrality"),
  make_option(c("-s", "--start-position"), type = "integer", default = 60000000,
              help = "Maximum start position for chr3p segments")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Function to extract GRanges from AneuFinder model
extract_granges <- function(file) {
  e <- new.env()
  load(file, envir = e)
  
  model <- NULL
  for (name in ls(e)) {
    if (is.list(e[[name]]) && "bins" %in% names(e[[name]])) {
      model <- e[[name]]
      break
    }
  }
  
  if (is.null(model)) {
    log_warn(paste0("No model object found in ", file))
    return(NULL)
  }
  
  return(model$bins)
}

# Main filtering function
filter_by_chr3p <- function() {
  log_info("Starting chr3p filtering")
  log_info(paste0("Input directory: ", opt$input))
  log_info(paste0("Output directory: ", opt$output))
  
  # Check if arm-level file exists
  if (!file.exists(opt$armlevel)) {
    log_error(paste0("Arm-level file not found: ", opt$armlevel))
    stop("Arm-level cytoband file required for filtering")
  }
  
  # Load arm-level information
  arms <- read.delim(opt$armlevel)
  arms$length <- arms$end - arms$start
  
  # Find AneuFinder model files
  model_dir <- file.path(opt$input, "MODELS", "method-edivisive")
  if (!dir.exists(model_dir)) {
    log_error(paste0("Model directory not found: ", model_dir))
    stop("AneuFinder models not found")
  }
  
  files <- list.files(model_dir, pattern = "\\.RData$", 
                     recursive = TRUE, full.names = TRUE)
  
  log_info(paste0("Found ", length(files), " model files"))
  
  # Extract sample names and GRanges
  sample_names <- gsub("\\.bam.*$", "", basename(files))
  gr_list <- setNames(lapply(files, extract_granges), sample_names)
  
  # Remove NULL entries (failed loads)
  gr_list <- gr_list[!sapply(gr_list, is.null)]
  
  if (length(gr_list) == 0) {
    log_error("No valid models could be loaded")
    stop("No valid AneuFinder models")
  }
  
  log_info(paste0("Successfully loaded ", length(gr_list), " models"))
  
  # Combine all GRanges
  gr_summary <- do.call(c, gr_list)
  
  # Calculate ploidy
  ploidy_raw <- sapply(gr_summary, function(gr) {
    sum(width(gr) * gr$copy.number) / sum(width(gr))
  })
  
  # Round ploidy
  ploidy_round <- sapply(ploidy_raw, function(p) {
    if (p > 3.3) return(4)
    if (p < 2 & p > 1.4) return(2)
    round(p)
  })
  
  # Calculate gain/loss
  gr_GainLoss <- mapply(function(gr, sample) {
    pl <- ploidy_round[sample]
    gr$ploidy <- pl
    gr$gainloss <- gr$copy.number - pl
    gr
  }, gr_summary, names(gr_summary), SIMPLIFY = FALSE)
  
  # Convert to data frame
  floris_GainLoss <- bind_rows(lapply(names(gr_GainLoss), function(sample) {
    gr <- gr_GainLoss[[sample]]
    data.frame(
      sample = sample,
      chrom = as.character(seqnames(gr)),
      start = start(gr),
      end = end(gr),
      strand = as.character(strand(gr)),
      counts = mcols(gr)$counts,
      mcounts = mcols(gr)$mcounts,
      pcounts = mcols(gr)$pcounts,
      state = mcols(gr)$state,
      copy_number = mcols(gr)$copy.number,
      gainloss = mcols(gr)$gainloss,
      ploidy = mcols(gr)$ploidy
    )
  }))
  
  # Split segments that span two arms
  log_info("Splitting segments spanning arm boundaries")
  rows_to_add <- list()
  
  for (i in 1:nrow(floris_GainLoss)) {
    test <- arms %>% filter(chrom == floris_GainLoss$chrom[i])
    
    if (nrow(test) >= 1 && test[1, 'end'] >= floris_GainLoss$start[i] && 
        test[1, 'end'] <= floris_GainLoss$end[i]) {
      # Segment spans centromere
      new_row <- floris_GainLoss[i, ]
      floris_GainLoss[i, 'end'] <- test[1, 'end']
      new_row['start'] <- test[1, 'end']
      rows_to_add[[length(rows_to_add) + 1]] <- new_row
    }
  }
  
  if (length(rows_to_add) > 0) {
    floris_GainLoss <- rbind(floris_GainLoss, do.call(rbind, rows_to_add))
  }
  
  # Add chromosome arm annotation
  log_info("Adding chromosome arm annotations")
  floris_GainLoss$chromarm <- NA
  
  for (i in 1:nrow(floris_GainLoss)) {
    test <- arms %>% filter(chrom == floris_GainLoss$chrom[i])
    
    if (nrow(test) >= 1) {
      if (floris_GainLoss$start[i] >= 0 && 
          floris_GainLoss$end[i] <= test[1, 'end']) {
        floris_GainLoss$chromarm[i] <- test[1, 'chromarms']
      } else if (nrow(test) >= 2 && 
                floris_GainLoss$start[i] >= test[2, 'start'] && 
                floris_GainLoss$end[i] <= test[2, 'end']) {
        floris_GainLoss$chromarm[i] <- test[2, 'chromarms']
      }
    }
  }
  
  # Merge segments with same chromarm, sample, and copy number
  log_info("Merging adjacent segments")
  floris_GainLoss_merged <- floris_GainLoss %>%
    arrange(sample, chromarm, start) %>%
    mutate(dup_ID = data.table::rleid(sample, chromarm, copy_number)) %>%
    group_by(sample, chromarm, copy_number, dup_ID) %>%
    summarise(
      chrom = first(chrom),
      start = first(start),
      end = last(end),
      ploidy = first(ploidy),
      gainloss = first(gainloss),
      state = first(state),
      counts = sum(counts),
      mcounts = sum(mcounts),
      pcounts = sum(pcounts),
      .groups = "drop"
    ) %>%
    select(-dup_ID)
  
  # Calculate segment lengths and proportions
  floris_GainLoss_merged <- floris_GainLoss_merged %>%
    mutate(seg_length = end - start)
  
  # Calculate proportion lengths
  for (i in 1:nrow(floris_GainLoss_merged)) {
    arm <- floris_GainLoss_merged$chromarm[i]
    if (!is.na(arm) && arm %in% arms$chromarms) {
      arm_length <- arms[arms$chromarms == arm, "length"]
      floris_GainLoss_merged[i, 'prop_length'] <- 
        floris_GainLoss_merged[i, "seg_length"] / arm_length
    }
  }
  
  # Filter based on chr3p neutrality
  log_info(paste0("Filtering based on ", opt$chromosome, " neutrality"))
  
  chr3p_neutral <- floris_GainLoss_merged %>%
    filter(prop_length > opt$`prop-threshold`) %>%
    filter(chromarm == opt$chromosome & 
           gainloss >= 0 & 
           start <= opt$`start-position`)
  
  # Create blacklist
  blacklist <- chr3p_neutral %>% 
    distinct(sample) %>%
    pull(sample)
  
  log_info(paste0("Identified ", length(blacklist), " samples for filtering"))
  
  # Create output directory
  dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
  
  # Filter and copy good samples
  filtered_samples <- setdiff(names(gr_list), blacklist)
  log_info(paste0("Keeping ", length(filtered_samples), " samples"))
  
  # Copy good model files to output
  output_model_dir <- file.path(opt$output, "MODELS", "method-edivisive")
  dir.create(output_model_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (sample in filtered_samples) {
    # Find original file
    orig_file <- files[grep(sample, files)]
    if (length(orig_file) > 0) {
      # Copy to output
      file.copy(orig_file[1], 
               file.path(output_model_dir, basename(orig_file[1])),
               overwrite = TRUE)
    }
  }
  
  # Write filtering summary
  summary_file <- file.path(opt$output, "filtering_summary.txt")
  summary_lines <- c(
    "CHR3P FILTERING SUMMARY",
    "=======================",
    paste0("Date: ", Sys.Date()),
    paste0("Input directory: ", opt$input),
    paste0("Total samples: ", length(gr_list)),
    paste0("Filtered samples: ", length(blacklist)),
    paste0("Retained samples: ", length(filtered_samples)),
    paste0("Filtering rate: ", 
           round(length(blacklist) / length(gr_list) * 100, 2), "%"),
    "",
    "FILTERING CRITERIA:",
    paste0("- Chromosome arm: ", opt$chromosome),
    paste0("- Proportion threshold: ", opt$`prop-threshold`),
    paste0("- Start position cutoff: ", format(opt$`start-position`, big.mark = ",")),
    paste0("- Gain/loss threshold: >= 0"),
    "",
    "FILTERED SAMPLES:",
    paste0("- ", blacklist)
  )
  
  writeLines(summary_lines, summary_file)
  log_info(paste0("Summary written to: ", summary_file))
  
  # Write filtered sample list
  filtered_list_file <- file.path(opt$output, "filtered_samples.txt")
  writeLines(blacklist, filtered_list_file)
  
  # Write retained sample list
  retained_list_file <- file.path(opt$output, "retained_samples.txt")
  writeLines(filtered_samples, retained_list_file)
  
  log_info("Filtering complete")
}

# Main execution
if (!interactive()) {
  tryCatch({
    filter_by_chr3p()
  }, error = function(e) {
    log_error(paste0("Fatal error: ", e$message))
    quit(status = 1)
  })
}
