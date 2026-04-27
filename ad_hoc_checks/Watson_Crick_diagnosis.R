# Node setup chunks
if (Sys.info()["nodename"] == "DDJ9Y7V9GC") {
  setwd("/Volumes/TracerX/")
} else {
  setwd("/camp/project/tracerX/")
}

library(tidyverse)
library(AneuFinder)
library(grid)
library(tools)

plate <- "plate17_4"

bam_dir <- paste0("working/VCAM1_GnT/DATA/384_well/plate17/", plate, "/bam/")
bam_files <- list.files(path = bam_dir, pattern = "*.bam$", full.names = TRUE, recursive = TRUE)

strand_plot_dir <- paste0("working/VCAM1_GnT/DATA/384_well/plate17/", plate, "/strandseq_plots")
dir.create(strand_plot_dir, recursive = TRUE, showWarnings = FALSE)

# containers for summary metrics
summary_bins_list <- list()
summary_sample_list <- list()

build_chr_guides <- function(df) {
  chr_pos_df <- df %>%
    group_by(seqnames) %>%
    summarise(
      start = min(cum_bin),
      end = max(cum_bin),
      .groups = "drop"
    ) %>%
    arrange(start)

  boundary_lines_df <- chr_pos_df[-nrow(chr_pos_df), , drop = FALSE]
  boundary_lines_df$xint <- boundary_lines_df$end + 0.5
  chr_pos_df$mid <- (chr_pos_df$start + chr_pos_df$end) / 2

  list(chr_pos_df = chr_pos_df, boundary_lines_df = boundary_lines_df)
}

make_prop_plot <- function(df_long, chr_pos_df, boundary_lines_df, sample_id) {
  ggplot(df_long, aes(x = cum_bin, y = prop, fill = strand)) +
    geom_col(width = 1) +
    geom_vline(
      data = boundary_lines_df,
      aes(xintercept = xint),
      linetype = "dashed",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      breaks = chr_pos_df$mid,
      labels = chr_pos_df$seqnames,
      expand = c(0, 0)
    ) +
    scale_fill_manual(values = c(plus = "#4477AA", minus = "#CC6677")) +
    labs(
      title = paste0("Genome wide strand proportions: ", sample_id),
      x = "Chromosome",
      y = "Proportion of reads",
      fill = "Strand"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.spacing.x = unit(0, "lines")
    )
}

make_strandseq_plot <- function(df, chr_pos_df, boundary_lines_df, sample_id, plate) {
  chr_lengths_df <- df %>%
    group_by(seqnames) %>%
    summarise(
      chr_end_mb = max(end) / 1e6,
      .groups = "drop"
    ) %>%
    arrange(factor(seqnames, levels = chr_pos_df$seqnames)) %>%
    mutate(chr_idx = seq_len(n()))

  strand_df <- df %>%
    transmute(
      seqnames = as.character(seqnames),
      start_mb = start / 1e6,
      end_mb = end / 1e6,
      watson = pcounts,
      crick = mcounts,
      total = total
    ) %>%
    left_join(chr_lengths_df %>% select(seqnames, chr_idx), by = "seqnames")

  max_bin_count <- max(c(strand_df$watson, strand_df$crick), na.rm = TRUE)
  max_bin_count <- ifelse(is.finite(max_bin_count) && max_bin_count > 0, max_bin_count, 1)

  center_gap <- 0.03
  max_half_width <- 0.42

  strand_df <- strand_df %>%
    mutate(
      watson_width = (watson / max_bin_count) * max_half_width,
      crick_width = (crick / max_bin_count) * max_half_width,
      watson_xmin = chr_idx - center_gap - watson_width,
      watson_xmax = chr_idx - center_gap,
      crick_xmin = chr_idx + center_gap,
      crick_xmax = chr_idx + center_gap + crick_width
    )

  annotation_text <- paste0(
    "Plate: ", plate,
    " | Sample: ", sample_id,
    " | Total reads: ", scales::comma(sum(df$total, na.rm = TRUE)),
    " | Median reads/bin: ", round(median(df$total, na.rm = TRUE), 1)
  )

  ggplot() +
    geom_rect(
      data = strand_df,
      aes(
        xmin = watson_xmin,
        xmax = watson_xmax,
        ymin = start_mb,
        ymax = end_mb
      ),
      fill = "#E6954A",
      colour = NA,
      alpha = 0.9
    ) +
    geom_rect(
      data = strand_df,
      aes(
        xmin = crick_xmin,
        xmax = crick_xmax,
        ymin = start_mb,
        ymax = end_mb
      ),
      fill = "#6E99A1",
      colour = NA,
      alpha = 0.9
    ) +
    geom_segment(
      data = chr_lengths_df,
      aes(
        x = chr_idx,
        xend = chr_idx,
        y = 0,
        yend = chr_end_mb
      ),
      linewidth = 0.35,
      colour = "black"
    ) +
    annotate(
      geom = "text",
      x = 1,
      y = max(chr_lengths_df$chr_end_mb, na.rm = TRUE) * 1.03,
      label = annotation_text,
      hjust = 0,
      size = 3.2
    ) +
    annotate(
      geom = "text",
      x = max(chr_lengths_df$chr_idx) / 2,
      y = -max(chr_lengths_df$chr_end_mb, na.rm = TRUE) * 0.05,
      label = "Watson | Crick",
      size = 3.8
    ) +
    scale_x_continuous(
      breaks = chr_lengths_df$chr_idx,
      labels = gsub("chr", "", chr_lengths_df$seqnames),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(-max(chr_lengths_df$chr_end_mb, na.rm = TRUE) * 0.08,
                 max(chr_lengths_df$chr_end_mb, na.rm = TRUE) * 1.06),
      labels = function(x) paste0(x, " Mb")
    ) +
    labs(
      title = paste0("Strand-seq Watson/Crick profile: ", sample_id),
      subtitle = "Chromosomes are vertical; Watson is plotted to the left and Crick to the right of each chromosome axis",
      x = "Chromosome",
      y = NULL
    ) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9),
      plot.title.position = "plot"
    )
}

## 1) Per BAM genome wide stacked plots (per page) ----------------------
pdf(paste0("working/VCAM1_GnT/DATA/384_well/plate17/", plate, "/strandseq_plots/p_mcounts_stacked.pdf"),
    width = 10, height = 6)

strand_plot_list <- list()

for (bam in bam_files) {
  message("Processing: ", bam)

  binned <- binReads(
    file = bam,
    binsize = 1e6,
    chromosomes = paste0("chr", c(1:22, "X", "Y")),
    assembly = "hg38",
    pairedEndReads = TRUE,
    remove.duplicate.reads = TRUE,
    min.mapq = 10
  )

  gr <- binned[[1]]
  df <- as.data.frame(gr)

  df$total <- df$pcounts + df$mcounts
  df <- df[df$total > 0, ]

  df <- df %>%
    arrange(seqnames, start)

  df$cum_bin <- seq_len(nrow(df))

  sample_id <- attr(gr, "ID")
  if (is.null(sample_id)) {
    sample_id <- file_path_sans_ext(basename(bam))
  }

  prop_plus_vec <- df$pcounts / df$total
  dev_vec <- abs(prop_plus_vec - 0.5)

  summary_bins_list[[sample_id]] <- tibble(
    sample_id = sample_id,
    prop_plus = prop_plus_vec,
    dev = dev_vec
  )

  P_tot <- sum(df$pcounts)
  M_tot <- sum(df$mcounts)
  n_tot <- P_tot + M_tot
  prop_plus_global <- P_tot / n_tot

  median_dev <- median(dev_vec, na.rm = TRUE)
  q90_dev <- as.numeric(quantile(dev_vec, 0.9, na.rm = TRUE))
  frac_dev_0.1 <- mean(dev_vec > 0.1, na.rm = TRUE)

  expected_sd <- sqrt(0.25 / n_tot)
  z_approx <- (prop_plus_global - 0.5) / expected_sd
  pval_approx <- 2 * pnorm(abs(z_approx), lower.tail = FALSE)

  n_vec <- df$total
  z_vec <- (prop_plus_vec - 0.5) / sqrt(0.25 / n_vec)

  chi2 <- sum(z_vec^2, na.rm = TRUE)
  df_chi <- length(z_vec) - 1
  phi <- chi2 / df_chi
  p_over <- pchisq(chi2, df = df_chi, lower.tail = FALSE)

  summary_sample_list[[sample_id]] <- tibble(
    sample_id = sample_id,
    prop_plus_global = prop_plus_global,
    n_tot = n_tot,
    median_dev = median_dev,
    q90_dev = q90_dev,
    frac_dev_0.1 = frac_dev_0.1,
    pval_approx = pval_approx,
    overdisp_phi = phi,
    overdisp_p = p_over
  )

  chr_guides <- build_chr_guides(df)
  chr_pos_df <- chr_guides$chr_pos_df
  boundary_lines_df <- chr_guides$boundary_lines_df

  df_long <- df %>%
    transmute(
      seqnames = as.character(seqnames),
      cum_bin,
      p = pcounts / total,
      m = mcounts / total
    ) %>%
    pivot_longer(
      cols = c(p, m),
      names_to = "strand",
      values_to = "prop"
    ) %>%
    mutate(strand = ifelse(strand == "p", "plus", "minus"))

  p <- make_prop_plot(df_long, chr_pos_df, boundary_lines_df, sample_id)
  print(p)

  strand_plot <- make_strandseq_plot(df, chr_pos_df, boundary_lines_df, sample_id, plate)
  strand_plot_list[[sample_id]] <- strand_plot

  ggsave(
    filename = file.path(strand_plot_dir, paste0(sample_id, "_strandseq_profile.pdf")),
    plot = strand_plot,
    width = 14,
    height = 5,
    units = "in"
  )
}

dev.off()

## 1b) Combined Watson/Crick strand-seq plots ---------------------------
pdf(paste0("working/VCAM1_GnT/DATA/384_well/plate17/", plate, "/strandseq_plots/watson_crick_strandseq_profiles.pdf"),
    width = 14, height = 5)

for (sample_id in names(strand_plot_list)) {
  print(strand_plot_list[[sample_id]])
}

dev.off()

## 2) Summary across samples (second PDF) -------------------------------
pdf(paste0("working/VCAM1_GnT/DATA/384_well/plate17/", plate, "/strandseq_plots/wat_crick_summary_diagnosis.pdf"),
    width = 14, height = 12)

summary_bins_df <- bind_rows(summary_bins_list)
summary_sample_df <- bind_rows(summary_sample_list)

summary_sample_df <- summary_sample_df %>%
  arrange(prop_plus_global)

summary_sample_df$sample_id <- factor(
  summary_sample_df$sample_id,
  levels = summary_sample_df$sample_id
)

summary_bins_df$sample_id <- factor(
  summary_bins_df$sample_id,
  levels = levels(summary_sample_df$sample_id)
)

summary_sample_long <- summary_sample_df %>%
  transmute(
    sample_id,
    plus = prop_plus_global,
    minus = 1 - prop_plus_global
  ) %>%
  pivot_longer(
    cols = c(plus, minus),
    names_to = "strand",
    values_to = "prop"
  )

p_global <- ggplot(summary_sample_long,
                   aes(y = sample_id, x = prop, fill = strand)) +
  geom_col() +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             colour = "red") +
  scale_x_continuous(limits = c(0, 1)) +
  labs(
    title = paste0("Global strand proportions per sample: ", plate),
    x = "Proportion of reads",
    y = "Sample"
  ) +
  theme_bw() +
  theme(legend.position = "top")

print(p_global)

p_violin <- ggplot(summary_bins_df,
                   aes(x = sample_id, y = prop_plus)) +
  geom_violin(scale = "width") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             colour = "red") +
  coord_flip() +
  labs(
    title = paste0("Per bin plus strand proportion: ", plate),
    x = "Sample",
    y = "Per bin plus proportion"
  ) +
  theme_bw()

print(p_violin)

p_bad <- ggplot(summary_sample_df,
                aes(x = sample_id, y = frac_dev_0.1)) +
  geom_col() +
  coord_flip() +
  labs(
    title = paste0("Fraction of bins with |plus - 0.5| > 0.1: ", plate),
    x = "Sample",
    y = "Fraction of bins"
  ) +
  theme_bw()

print(p_bad)

summary_sample_df <- summary_sample_df %>%
  mutate(
    overdisp_signif = case_when(
      overdisp_p < 1e-6 ~ "***",
      overdisp_p < 1e-3 ~ "**",
      overdisp_p < 0.01 ~ "*",
      TRUE ~ ""
    )
  )

max_phi <- max(summary_sample_df$overdisp_phi, na.rm = TRUE)

p_overdisp <- ggplot(summary_sample_df,
                     aes(x = sample_id, y = overdisp_phi)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  geom_text(aes(
    label = overdisp_signif,
    y = overdisp_phi + 0.05 * max_phi
  ),
  hjust = 0) +
  expand_limits(y = max_phi * 1.15) +
  coord_flip() +
  labs(
    title = paste0("Overdispersion of strand proportions: ", plate),
    x = "Sample",
    y = "Overdispersion factor (phi)"
  ) +
  theme_bw()

print(p_overdisp)

dev.off()
