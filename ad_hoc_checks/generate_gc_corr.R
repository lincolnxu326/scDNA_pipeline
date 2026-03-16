#!/usr/bin/env Rscript

library(AneuFinder)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)

plate_dir <- "plate7_pipe_test" #EDIT here

VARIABLE_WIDTH_REFERENCE <- "/nemo/project/proj-tracerX/working/VCAM1_GnT/DAQI/source/Aneufinder/mappability/normal_ref_sorted.bam.bed"
BLACKLIST <-paste0("/nemo/project/proj-tracerX/working/VCAM1_GnT/DATA/",plate_dir,"/mappability/blacklist.bed.gz")   # set to the same file you use in the real run, or leave NULL
CHROMS <- paste0("chr", c(1:22))

dir.create(paste0("/nemo/project/proj-tracerX/working/VCAM1_GnT/DATA/",plate_dir,"/GC"),recursive = T)


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
    gc_vals[idx] <- rowSums(freqs[, c("G","C"), drop = FALSE])
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
      gc_vals[idx] <- rowSums(freqs[, c("G","C"), drop = FALSE])
    }
    mcols(gr)$GC <- gc_vals
    ref_binned_[[i]] <- gr
  }
}

saveRDS(
  ref_binned_,
  paste0("/nemo/project/proj-tracerX/working/VCAM1_GnT/DATA/",plate_dir,"/GC/hg38_variable_bins_with_GC.rds"),version=2
)