#!/usr/bin/env Rscript

if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

options(download.file.method = "wget")

suppressPackageStartupMessages({
  aneufinder_available <- requireNamespace("AneuFinder", quietly = TRUE)
  if (aneufinder_available) {
    library(AneuFinder)
  } else {
    library(GenomicRanges)
    library(Rsamtools)
  }
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Input mappability BAM file (REQUIRED)"),
  make_option(c("-b", "--binsize"), type = "integer", default = 100000,
              help = "Bin size for mappability (default: 100kb)"),
  make_option(c("-L", "--lower-quantile"), type = "numeric", default = 0.12,
              help = "Lower quantile for blacklist (default: 0.12)"),
  make_option(c("-U", "--upper-quantile"), type = "numeric", default = 0.95,
              help = "Upper quantile for blacklist (default: 0.95)"),
  make_option(c("-c", "--chromosomes"), type = "character",
              default = paste0("chr", c(1:22), collapse = ","),
              help = "Chromosomes to analyze (comma-separated, default: chr1-22)"),
  make_option(c("--output-bed"), type = "character",
              help = "Output BED file for blacklist (REQUIRED)"),
  make_option(c("--output-plot"), type = "character",
              help = "Output PDF for diagnosis plot"),
  make_option(c("--output-stats"), type = "character",
              help = "Output text file with statistics")
)

opt_parser <- OptionParser(option_list = option_list,
                           description = "Generate blacklist from mappability BAM file")
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  stop("Input BAM file is required. Use --input to specify the mappability BAM file.")
}
if (is.null(opt$`output-bed`)) {
  stop("Output BED file is required. Use --output-bed to specify the output blacklist.")
}

generate_blacklist_from_bam <- function() {
  cat("========================================\n")
  cat("Generating blacklist from mappability BAM\n")
  cat("========================================\n")
  cat(paste0("Input BAM: ", opt$input, "\n"))
  
  if (!file.exists(opt$input)) {
    stop(paste0("Input BAM file not found: ", opt$input))
  }
  
  chromosomes <- unlist(strsplit(opt$chromosomes, ","))
  cat(paste0("Chromosomes: ", paste(chromosomes, collapse = ", "), "\n"))
  cat(paste0("Binning reads with bin size: ", opt$binsize, " bp\n"))
  
  if (aneufinder_available) {
    # back to positional arg, because your version does not know filename=
    bins <- binReads(
      opt$input,
      assembly = "hg38",
      binsizes = opt$binsize,
      chromosomes = chromosomes
    )[[1]]
    counts <- bins$counts
  } else {
    cat("Using fallback method (AneuFinder not available)\n")
    bam <- BamFile(opt$input)
    header <- scanBamHeader(bam)
    chrom_sizes <- header[[1]]$targets
    chrom_sizes <- chrom_sizes[names(chrom_sizes) %in% chromosomes]
    if (length(chrom_sizes) == 0) {
      stop("No matching chromosomes found in BAM file")
    }
    all_bins <- GRanges()
    for (chr in names(chrom_sizes)) {
      starts_ <- seq(1, chrom_sizes[chr], by = opt$binsize)
      chr_bins <- GRanges(
        seqnames = chr,
        ranges = IRanges(
          start = starts_,
          width = pmin(opt$binsize, chrom_sizes[chr] - starts_ + 1)
        )
      )
      all_bins <- c(all_bins, chr_bins)
    }
    cat("Counting reads in bins...\n")
    param <- ScanBamParam(which = all_bins)
    bin_counts <- countBam(bam, param = param)
    bins <- all_bins
    mcols(bins)$counts <- bin_counts$records
    counts <- bin_counts$records
  }
  
  cat(paste0("Total bins: ", length(bins), "\n"))
  
  lcutoff <- quantile(counts, opt$`lower-quantile`)
  ucutoff <- quantile(counts, opt$`upper-quantile`)
  cat(paste0("Lower cutoff (", opt$`lower-quantile`, " quantile): ", round(lcutoff, 2), "\n"))
  cat(paste0("Upper cutoff (", opt$`upper-quantile`, " quantile): ", round(ucutoff, 2), "\n"))
  
  if (!is.null(opt$`output-plot`)) {
    cat(paste0("Creating diagnosis plot: ", opt$`output-plot`, "\n"))
    bins_df <- as.data.frame(bins)
    bins_df$counts <- counts
    bins_df$chr <- factor(as.character(seqnames(bins)), levels = chromosomes)
    
    p <- ggplot(bins_df, aes(x = start/1e6, y = counts)) +
      geom_point(size = 0.5, alpha = 0.5) +
      facet_wrap(~ chr, scales = "free_x", ncol = 4) +
      geom_hline(yintercept = lcutoff, color = "red", linetype = "dashed", size = 0.5, alpha = 0.7) +
      geom_hline(yintercept = ucutoff, color = "red", linetype = "dashed", size = 0.5, alpha = 0.7) +
      theme_minimal() +
      labs(
        title = "Mappability Analysis for Blacklist Generation",
        subtitle = paste0("Bin size: ", opt$binsize/1000, "kb | Lower cutoff: ", round(lcutoff, 1),
                          " | Upper cutoff: ", round(ucutoff, 1)),
        x = "Position (Mb)",
        y = "Read count per bin"
      )
    
    p_hist <- ggplot(bins_df, aes(x = counts)) +
      geom_histogram(bins = 50, fill = "gray40", alpha = 0.7) +
      geom_vline(xintercept = lcutoff, color = "red", linetype = "dashed") +
      geom_vline(xintercept = ucutoff, color = "red", linetype = "dashed") +
      theme_minimal() +
      labs(
        title = "Distribution of bin counts",
        subtitle = paste0("Total bins: ", length(bins)),
        x = "Counts",
        y = "Frequency"
      )
    
    pdf(opt$`output-plot`, width = 12, height = 10)
    print(p)
    print(p_hist)
    dev.off()
  }
  
  blacklist <- bins[counts <= lcutoff | counts >= ucutoff]
  blacklist <- reduce(blacklist)
  cat(paste0("Number of blacklisted regions: ", length(blacklist), "\n"))
  cat(paste0("Total blacklisted bases: ", format(sum(width(blacklist)), big.mark = ","), " bp\n"))
  
  cat(paste0("Writing blacklist to: ", opt$`output-bed`, "\n"))
  path_plain_ <- sub("\\.gz$", "", opt$`output-bed`)
  
  if (aneufinder_available) {
    # AneuFinder writes path_plain_.bed.gz
    exportGRanges(
      blacklist,
      filename = path_plain_,
      header = FALSE,
      chromosome.format = "NCBI"
    )
    exported_file_ <- paste0(path_plain_, ".bed.gz")
    
    if (grepl("\\.gz$", opt$`output-bed`)) {
      if (file.exists(exported_file_)) {
        file.rename(exported_file_, opt$`output-bed`)
      } else {
        stop("Expected AneuFinder output not found: ", exported_file_)
      }
    } else {
      if (file.exists(exported_file_)) {
        system(paste0("gunzip -f ", exported_file_))
      } else {
        stop("Expected AneuFinder output not found: ", exported_file_)
      }
    }
    
  } else {
    # fallback writer
    blacklist_df <- data.frame(
      chrom = as.character(seqnames(blacklist)),
      start = start(blacklist) - 1,
      end   = end(blacklist)
    )
    write.table(
      blacklist_df,
      file = path_plain_,
      col.names = FALSE,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    if (grepl("\\.gz$", opt$`output-bed`)) {
      system(paste0("gzip -f ", path_plain_))
    }
  }
  
  if (!is.null(opt$`output-stats`)) {
    cat(paste0("Writing statistics to: ", opt$`output-stats`, "\n"))
    stats_text <- c(
      "BLACKLIST GENERATION STATISTICS",
      "================================",
      paste0("Date: ", Sys.Date()),
      paste0("Input BAM: ", opt$input),
      paste0("Bin size: ", format(opt$binsize, big.mark = ","), " bp"),
      paste0("Chromosomes analyzed: ", paste(chromosomes, collapse = ", ")),
      "",
      "CUTOFF VALUES:",
      paste0("Lower quantile: ", opt$`lower-quantile`, " (cutoff: ", round(lcutoff, 2), ")"),
      paste0("Upper quantile: ", opt$`upper-quantile`, " (cutoff: ", round(ucutoff, 2), ")")
    )
    writeLines(stats_text, opt$`output-stats`)
  }
  
  cat("\nBlacklist generation complete!\n")
}

main <- function() {
  tryCatch({
    generate_blacklist_from_bam()
  }, error = function(e) {
    cat("\nERROR: ", e$message, "\n", sep = "")
    quit(status = 1)
  })
}

if (!interactive()) {
  main()
}
