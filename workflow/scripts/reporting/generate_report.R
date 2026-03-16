#!/usr/bin/env Rscript
#' Generate HTML report for single-cell sequencing pipeline
#' 
#' Creates comprehensive HTML report with QC metrics and visualizations

# Activate renv for reproducible R environment
if (file.exists("renv/activate.R")) {
  source("renv/activate.R")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(jsonlite)
  library(optparse)
  library(knitr)
  library(rmarkdown)
})

# Parse arguments
option_list <- list(
  make_option(c("-j", "--qc-json"), type = "character",
              help = "QC JSON file from compile_qc.py"),
  make_option(c("-t", "--qc-tsv"), type = "character",
              help = "QC TSV summary table"),
  make_option(c("-a", "--aneufinder-dir"), type = "character",
              help = "AneuFinder output directory"),
  make_option(c("-o", "--output"), type = "character",
              help = "Output HTML report path"),
  make_option(c("-p", "--plate"), type = "character",
              help = "Plate name")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Function to create report
generate_report <- function() {
  # Load QC data
  qc_data <- fromJSON(opt$`qc-json`)
  qc_summary <- read.table(opt$`qc-tsv`, header = TRUE, sep = "\t")
  
  # Create temporary Rmd file
  rmd_content <- '
---
title: "Single-Cell DNA Sequencing Report"
subtitle: "Plate: `r plate_name`"
date: "`r format(Sys.Date(), "%B %d, %Y")`"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float: true
    toc_depth: 3
    number_sections: false
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, 
                      fig.width = 10, fig.height = 6)

library(tidyverse)
library(ggplot2)
library(plotly)
library(DT)

# Set theme
theme_set(theme_minimal())
```

# Overview {.tabset}

## Summary Statistics

```{r summary}
# Create summary cards
summary_stats <- data.frame(
  Metric = c("Total Reads", "Assigned Reads", "Unique Molecules", 
             "Wells Processed", "Average Duplication Rate", "Average Mapping Rate"),
  Value = c(
    format(qc_data$summary$total_reads, big.mark = ","),
    format(qc_data$summary$assigned_reads, big.mark = ","),
    format(qc_data$summary$total_unique, big.mark = ","),
    qc_data$summary$wells_processed,
    paste0(round(qc_data$summary$avg_duplication_rate, 1), "%"),
    paste0(round(qc_data$summary$avg_mapping_rate, 1), "%")
  )
)

DT::datatable(summary_stats, 
              options = list(dom = "t", paging = FALSE),
              rownames = FALSE)
```

## Quality Metrics

```{r metrics-plot}
# Create quality metrics plot
metrics_df <- qc_summary %>%
  select(well_id, duplication_rate, mapping_rate, adapter_contamination) %>%
  pivot_longer(cols = -well_id, names_to = "metric", values_to = "value")

p <- ggplot(metrics_df, aes(x = well_id, y = value, fill = metric)) +
  geom_col() +
  facet_wrap(~ metric, scales = "free_y", ncol = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  labs(title = "Quality Metrics by Well",
       x = "Well ID",
       y = "Value (%)") +
  scale_fill_viridis_d()

ggplotly(p)
```

# Demultiplexing Results

```{r demux}
# Demultiplexing statistics
demux_stats <- data.frame(
  well_id = names(qc_data$demux$per_well),
  reads = unlist(qc_data$demux$per_well)
) %>%
  arrange(desc(reads))

p_demux <- ggplot(demux_stats, aes(x = reorder(well_id, reads), y = reads)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Read Distribution Across Wells",
       x = "Well ID",
       y = "Number of Reads") +
  theme(axis.text.y = element_text(size = 6))

ggplotly(p_demux)
```

## Assignment Rate

```{r assignment}
assignment_data <- data.frame(
  Category = c("Assigned", "Unassigned"),
  Reads = c(qc_data$summary$assigned_reads, 
           qc_data$summary$total_reads - qc_data$summary$assigned_reads)
)

p_assign <- ggplot(assignment_data, aes(x = "", y = Reads, fill = Category)) +
  geom_col(width = 1) +
  coord_polar("y") +
  labs(title = "Barcode Assignment Rate") +
  scale_fill_manual(values = c("Assigned" = "#2ecc71", "Unassigned" = "#e74c3c"))

ggplotly(p_assign)
```

# Deduplication Analysis

```{r dedup}
# Deduplication rates by well
dedup_df <- qc_summary %>%
  select(well_id, total_reads, unique_reads, duplication_rate) %>%
  arrange(desc(duplication_rate))

p_dedup <- ggplot(dedup_df, aes(x = reorder(well_id, -duplication_rate), 
                                y = duplication_rate)) +
  geom_col(aes(fill = duplication_rate)) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
  scale_fill_gradient2(low = "green", mid = "yellow", high = "red", midpoint = 50) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(title = "PCR Duplication Rates by Well",
       subtitle = "Red line indicates 50% threshold",
       x = "Well ID",
       y = "Duplication Rate (%)")

ggplotly(p_dedup)
```

# Alignment Statistics

```{r alignment}
# Alignment rates
align_df <- qc_summary %>%
  filter(!is.na(mapping_rate)) %>%
  select(well_id, mapping_rate) %>%
  arrange(desc(mapping_rate))

p_align <- ggplot(align_df, aes(x = reorder(well_id, -mapping_rate), 
                                y = mapping_rate)) +
  geom_col(aes(fill = mapping_rate)) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "red") +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "green", midpoint = 70) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(title = "Alignment Rates by Well",
       subtitle = "Red line indicates 70% threshold",
       x = "Well ID",
       y = "Mapping Rate (%)")

ggplotly(p_align)
```

# Copy Number Analysis

```{r cnv, eval=file.exists(opt$aneufinder_dir)}
# Check if AneuFinder results exist
aneufinder_available <- FALSE
if (!is.null(opt$aneufinder_dir) && dir.exists(opt$aneufinder_dir)) {
  model_files <- list.files(file.path(opt$aneufinder_dir, "MODELS"),
                           pattern = "\\\\.RData$", 
                           recursive = TRUE, full.names = TRUE)
  if (length(model_files) > 0) {
    aneufinder_available <- TRUE
    cat(paste0("Found ", length(model_files), " AneuFinder models\\n"))
  }
}

if (!aneufinder_available) {
  cat("AneuFinder results not available\\n")
}
```

**Note:** GC correction is currently disabled due to BSgenome.Hsapiens.UCSC.hg38 not being available in the environment.

# Quality Control Summary

## Pass/Fail Summary

```{r qc-summary}
# Apply QC thresholds
qc_summary$pass_duplication <- qc_summary$duplication_rate < 50
qc_summary$pass_mapping <- qc_summary$mapping_rate > 70
qc_summary$pass_adapter <- qc_summary$adapter_contamination < 10
qc_summary$pass_all <- qc_summary$pass_duplication & 
                       qc_summary$pass_mapping & 
                       qc_summary$pass_adapter

# Summary table
pass_summary <- data.frame(
  Criterion = c("Duplication Rate < 50%", "Mapping Rate > 70%", 
                "Adapter Contamination < 10%", "Pass All Criteria"),
  Pass = c(sum(qc_summary$pass_duplication, na.rm = TRUE),
          sum(qc_summary$pass_mapping, na.rm = TRUE),
          sum(qc_summary$pass_adapter, na.rm = TRUE),
          sum(qc_summary$pass_all, na.rm = TRUE)),
  Fail = c(sum(!qc_summary$pass_duplication, na.rm = TRUE),
          sum(!qc_summary$pass_mapping, na.rm = TRUE),
          sum(!qc_summary$pass_adapter, na.rm = TRUE),
          sum(!qc_summary$pass_all, na.rm = TRUE))
)

DT::datatable(pass_summary, 
              options = list(dom = "t", paging = FALSE),
              rownames = FALSE) %>%
  DT::formatStyle("Pass", backgroundColor = styleInterval(0, c("white", "#d4edda"))) %>%
  DT::formatStyle("Fail", backgroundColor = styleInterval(0, c("white", "#f8d7da")))
```

## Detailed Well Status

```{r well-status}
# Create detailed status table
well_status <- qc_summary %>%
  select(well_id, total_reads, duplication_rate, mapping_rate, 
         adapter_contamination, pass_all) %>%
  mutate(Status = ifelse(pass_all, "PASS", "FAIL"))

DT::datatable(well_status,
              options = list(pageLength = 10),
              rownames = FALSE) %>%
  DT::formatStyle("Status",
                 backgroundColor = styleEqual(c("PASS", "FAIL"),
                                             c("#d4edda", "#f8d7da")))
```

# Pipeline Information

```{r pipeline-info}
pipeline_info <- data.frame(
  Parameter = c("Pipeline Version", "Plate", "Date Processed", 
                "Reference Genome", "Bin Size", "Segmentation Method"),
  Value = c("2.0", plate_name, format(Sys.Date(), "%Y-%m-%d"),
           "GRCh38", "1 Mb", "edivisive")
)

DT::datatable(pipeline_info, 
              options = list(dom = "t", paging = FALSE),
              rownames = FALSE)
```

---

*Report generated automatically by scDNA sequencing pipeline*
'
  
  # Write Rmd file
  rmd_file <- tempfile(fileext = ".Rmd")
  writeLines(rmd_content, rmd_file)
  
  # Render report
  rmarkdown::render(
    input = rmd_file,
    output_file = opt$output,
    params = list(
      plate_name = opt$plate
    ),
    envir = new.env()
  )
  
  # Clean up
  unlink(rmd_file)
  
  cat(paste0("Report generated: ", opt$output, "\n"))
}

# Main execution
if (!interactive()) {
  tryCatch({
    generate_report()
  }, error = function(e) {
    cat(paste0("Error generating report: ", e$message, "\n"))
    quit(status = 1)
  })
}
