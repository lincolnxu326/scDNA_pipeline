library(QDNAseq)
library(QDNAseq.hg38)

# activate X11 display in Rstudio Server
#options(bitmapType='cairo')

setwd("/Volumes/TracerX/working/VCAM1_GnT")

bam_dir <- "/Volumes/TracerX/working/VCAM1_GnT/DATA/plate7_pipe_test/bam/"

dir.create("/Volumes/TracerX/working/VCAM1_GnT/DATA/plate10/aneufinder/per_bin_plots/",recursive = T)
output_dir <- "/Volumes/TracerX/working/VCAM1_GnT/DATA/plate10/aneufinder/per_bin_plots/"

# using simple mapping bams
# bam_dir <- "/nemo/project/tracerX/working/VCAM1_GnT/DAQI/output/plate7/simple_mapping/bam/"
# output_dir <- "/nemo/project/tracerX/working/VCAM1_GnT/DAQI/working/plate7/cnv_raw_1000kbb/"


bam_files <- list.files(path=bam_dir, pattern="*.bam$", full.names=TRUE, recursive=TRUE)


binSize <- 1000  # Use 100 kb bins (adjust based on data)
bins <- getBinAnnotations(binSize, genome = "hg38")

# readCounts <- binReadCounts(bins, bamfiles=paste0(bam_dir,"/W93.bam"))
# readCountsFiltered <- applyFilters(readCounts)
# readCountsFiltered <- estimateCorrection(readCountsFiltered)
# readCountsFiltered <- correctBins(readCountsFiltered)
# copyNumbers <- segmentBins(readCountsFiltered)
# p <- plot(readCountsFiltered, main="Copy Number Profile")

# do it in a lopp 
for (bam in bam_files) {
  # Extract sample name from BAM file path
  sample_name <- tools::file_path_sans_ext(basename(bam))
  # Bin read counts
  readCounts <- binReadCounts(bins, bamfiles=bam)
  # Apply filters and normalize
  readCountsFiltered <- applyFilters(readCounts)
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  readCountsFiltered <- correctBins(readCountsFiltered)
  # Segment copy number
  copyNumbers <- segmentBins(readCountsFiltered)
  
  # Save CNV plot
  png(filename=paste0(output_dir, sample_name, "_CNV.png"), width=2400, height=1000)
  plot(copyNumbers, main=paste("Copy Number Profile -", sample_name))
  dev.off()
}



