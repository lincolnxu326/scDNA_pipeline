#!/usr/bin/env Rscript

library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)

plate_dir <- "plate17/plate17_1" # EDIT here

VARIABLE_WIDTH_REFERENCE <- "/nemo/project/proj-tracerX/working/VCAM1_GnT/DAQI/source/Aneufinder/mappability/normal_ref_sorted.bam.bed"
BLACKLIST <- "/nemo/project/proj-tracerX/working/VCAM1_GnT/TOOLS/scDNA_pipeline/resources/reference/mappability/blacklist.bed.gz"
CHROMS <- paste0("chr", c(1:22))
OUTPUT_RDS <- "/nemo/project/proj-tracerX/working/VCAM1_GnT/TOOLS/scDNA_pipeline/resources/reference/hg38_binsize1000000_variable_bins_with_GC.rds"

dir.create(dirname(OUTPUT_RDS), recursive = TRUE, showWarnings = FALSE)

ref_binned_ <- binReads(
  file = VARIABLE_WIDTH_REFERENCE,
  chromosomes = CHROMS,
  binsizes = 1e6,
  variable.width.reference = VARIABLE_WIDTH_REFERENCE,
  blacklist = BLACKLIST,
  assembly = "hg38",
  save.as.RData = FALSE,
  use.bamsignals = FALSE
)

if (!(inherits(ref_binned_, "GRanges") || inherits(ref_binned_, "GRangesList"))) {
  ref_binned_ <- ref_binned_[[1]]
}

if (inherits(ref_binned_, "GRanges")) {
  gr <- ref_binned_
  gc_vals <- numeric(length(gr))
  for (chr in unique(as.character(seqnames(gr)))) {
    idx <- which(seqnames(gr) == chr)
    views <- Biostrings::Views(BSgenome.Hsapiens.UCSC.hg38[[chr]], ranges(gr)[idx])
    freqs <- Biostrings::alphabetFrequency(views, as.prob = TRUE, baseOnly = TRUE)
    gc_vals[idx] <- rowSums(freqs[, c("G", "C"), drop = FALSE])
  }
  mcols(gr)$GC <- gc_vals
  ref_binned_ <- GenomicRanges::GRangesList("0" = gr)
} else if (inherits(ref_binned_, "GRangesList")) {
  for (i in seq_along(ref_binned_)) {
    gr <- ref_binned_[[i]]
    gc_vals <- numeric(length(gr))
    for (chr in unique(as.character(seqnames(gr)))) {
      idx <- which(seqnames(gr) == chr)
      views <- Biostrings::Views(BSgenome.Hsapiens.UCSC.hg38[[chr]], ranges(gr)[idx])
      freqs <- Biostrings::alphabetFrequency(views, as.prob = TRUE, baseOnly = TRUE)
      gc_vals[idx] <- rowSums(freqs[, c("G", "C"), drop = FALSE])
    }
    mcols(gr)$GC <- gc_vals
    ref_binned_[[i]] <- gr
  }
}

saveRDS(ref_binned_, OUTPUT_RDS, version = 2)
