# Node setup chunks
if(Sys.info()["nodename"]=="DDJ9Y7V9GC"){
  setwd("/Volumes/TracerX/")
}else{
  setwd("/camp/project/tracerX/")
}

library(tidyverse)
library(AneuFinder)
library(grid)

plate <- "plate7_pipe_test"

bam_dir <- paste0("working/VCAM1_GnT/DATA/",plate,"/bam/")
bam_files <- list.files(path = bam_dir, pattern = "*.bam$", full.names = TRUE, recursive = TRUE)

# containers for summary metrics
summary_bins_list   <- list()
summary_sample_list <- list()



## 1) Per BAM genome wide stacked plots (per page) ----------------------
pdf(paste0("working/VCAM1_GnT/DATA/",plate,"/aneufinder/p_mcounts_stacked.pdf"),
    width = 10, height = 6)

for (bam in bam_files) {
  
  message("Processing: ", bam)
  
  binned <- binReads(
    file                   = bam,
    binsize                = 1e6,                     # tune if needed
    chromosomes            = paste0("chr", c(1:22, "X", "Y")),
    assembly               = "hg38",
    pairedEndReads         = TRUE,                    # FALSE if SE
    remove.duplicate.reads = TRUE,
    min.mapq               = 10
  )
  
  gr <- binned[[1]]
  df <- as.data.frame(gr)
  
  df$total <- df$pcounts + df$mcounts
  df <- df[df$total > 0, ]
  
  # order bins genome wide and make a cumulative bin index
  df <- df %>%
    arrange(seqnames, start)
  
  df$cum_bin <- seq_len(nrow(df))
  
  # sample ID
  sample_id <- attr(gr, "ID")
  if (is.null(sample_id)) sample_id <- file_path_sans_ext(basename(bam))
  
  # ---------- per-bin summary stats and QC for this sample -------------
  # Per bin plus strand proportion and deviation from 0.5
  prop_plus_vec <- df$pcounts / df$total
  dev_vec       <- abs(prop_plus_vec - 0.5)
  
  # Keep per bin values (for violin / boxplot later)
  summary_bins_list[[sample_id]] <- tibble(
    sample_id = sample_id,
    prop_plus = prop_plus_vec,
    dev       = dev_vec
  )
  
  # Global (summed) strand counts and simple global deviation
  P_tot <- sum(df$pcounts)
  M_tot <- sum(df$mcounts)
  n_tot <- P_tot + M_tot
  prop_plus_global <- P_tot / n_tot
  
  median_dev   <- median(dev_vec, na.rm = TRUE)
  q90_dev      <- as.numeric(quantile(dev_vec, 0.9, na.rm = TRUE))
  frac_dev_0.1 <- mean(dev_vec > 0.1, na.rm = TRUE)
  
  # Approximate binomial p value for global plus proportion
  expected_sd <- sqrt(0.25 / n_tot)
  z_approx    <- (prop_plus_global - 0.5) / expected_sd
  pval_approx <- 2 * pnorm(abs(z_approx), lower.tail = FALSE)
  
  # ---- Overdispersion QC (binomial model, p = 0.5 per bin) -----------
  # Expected variance per bin: Var(prop_plus) ≈ 0.25 / n
  n_vec <- df$total
  
  # Standardised residuals; should be N(0,1) if only binomial noise
  z_vec <- (prop_plus_vec - 0.5) / sqrt(0.25 / n_vec)
  # same as: z_vec <- 2 * (prop_plus_vec - 0.5) * sqrt(n_vec)
  
  # Chi squared statistic across bins
  chi2   <- sum(z_vec^2, na.rm = TRUE)
  df_chi <- length(z_vec) - 1
  
  # Overdispersion factor: phi ~ 1 if variance ~ expectation
  phi    <- chi2 / df_chi
  
  # P value for overdispersion (large chi2 -> small p)
  p_over <- pchisq(chi2, df = df_chi, lower.tail = FALSE)
  
  # Per sample summary row (store everything here)
  summary_sample_list[[sample_id]] <- tibble(
    sample_id        = sample_id,
    prop_plus_global = prop_plus_global,   # global plus fraction
    n_tot            = n_tot,              # total reads across genome
    median_dev       = median_dev,         # median |p - 0.5| across bins
    q90_dev          = q90_dev,            # 90th percentile deviation
    frac_dev_0.1     = frac_dev_0.1,       # fraction of bins outside [0.4, 0.6]
    pval_approx      = pval_approx,        # global deviation p value
    overdisp_phi     = phi,                # overdispersion factor
    overdisp_p       = p_over              # overdispersion p value
  )
  # ---------------------------------------------------------------------
  
  # chromosome boundaries and midpoints for ticks / vlines
  chr_pos_df <- df %>%
    group_by(seqnames) %>%
    summarise(
      start = min(cum_bin),
      end   = max(cum_bin),
      .groups = "drop"
    ) %>%
    arrange(start)
  
  boundary_lines_df <- chr_pos_df[-nrow(chr_pos_df), , drop = FALSE]
  boundary_lines_df$xint <- boundary_lines_df$end + 0.5
  
  chr_pos_df$mid <- (chr_pos_df$start + chr_pos_df$end) / 2
  
  # long format, plus/minus proportions (stacked to 1)
  df_long <- df %>%
    transmute(
      seqnames = as.character(seqnames),
      cum_bin,
      p = pcounts / total,
      m = mcounts / total
    ) %>%
    pivot_longer(
      cols      = c(p, m),
      names_to  = "strand",
      values_to = "prop"
    ) %>%
    mutate(
      strand = ifelse(strand == "p", "plus", "minus")
    )
  
  p <- ggplot(df_long, aes(x = cum_bin, y = prop, fill = strand)) +
    geom_col(width = 1) +
    geom_vline(
      data        = boundary_lines_df,
      aes(xintercept = xint),
      linetype    = "dashed",
      linewidth   = 0.3,
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      breaks = chr_pos_df$mid,
      labels = chr_pos_df$seqnames,
      expand = c(0, 0)
    ) +
    labs(
      title = paste0("Genome wide strand proportions: ", sample_id),
      x     = "Chromosome",
      y     = "Proportion of reads"
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.spacing.x = unit(0, "lines")
    )
  
  print(p)   # one page per BAM
}

dev.off()

## 2) Summary across samples (second PDF) -------------------------------

pdf(paste0("working/VCAM1_GnT/DATA/", plate,"/aneufinder/wat_crick_summary_diagnosis.pdf"),width = 14, height = 12)

summary_bins_df   <- bind_rows(summary_bins_list)
summary_sample_df <- bind_rows(summary_sample_list)

# order samples by global plus proportion
summary_sample_df <- summary_sample_df %>%
  arrange(prop_plus_global)

summary_sample_df$sample_id <- factor(summary_sample_df$sample_id,
                                      levels = summary_sample_df$sample_id)

summary_bins_df$sample_id <- factor(summary_bins_df$sample_id,
                                    levels = levels(summary_sample_df$sample_id))

# 1) global mean p/m per sample (stacked, with 50% line)
summary_sample_long <- summary_sample_df %>%
  transmute(
    sample_id,
    plus  = prop_plus_global,
    minus = 1 - prop_plus_global
  ) %>%
  pivot_longer(
    cols      = c(plus, minus),
    names_to  = "strand",
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
    x     = "Proportion of reads",
    y     = "Sample"
  ) +
  theme_bw() +
  theme(legend.position = "top")

print(p_global)

# 2) distribution of per bin plus proportion per sample
p_violin <- ggplot(summary_bins_df,
                   aes(x = sample_id, y = prop_plus)) +
  geom_violin(scale = "width") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             colour = "red") +
  coord_flip() +
  labs(
    title = paste0("Per bin plus strand proportion: ", plate),
    x     = "Sample",
    y     = "Per bin plus proportion"
  ) +
  theme_bw()

print(p_violin)

# 3) fraction of strongly unbalanced bins per sample (|p - 0.5| > 0.1)
p_bad <- ggplot(summary_sample_df,
                aes(x = sample_id, y = frac_dev_0.1)) +
  geom_col() +
  coord_flip() +
  labs(
    title = paste0("Fraction of bins with |plus - 0.5| > 0.1: ", plate),
    x     = "Sample",
    y     = "Fraction of bins"
  ) +
  theme_bw()

print(p_bad)

# 4) Overdispersion factor per sample with significance stars
#    Define a simple star code based on overdisp_p
summary_sample_df <- summary_sample_df %>%
  mutate(
    overdisp_signif = case_when(
      overdisp_p < 1e-6  ~ "***",
      overdisp_p < 1e-3  ~ "**",
      overdisp_p < 0.01  ~ "*",
      TRUE               ~ ""
    )
  )

max_phi <- max(summary_sample_df$overdisp_phi, na.rm = TRUE)

p_overdisp <- ggplot(summary_sample_df,
                     aes(x = sample_id, y = overdisp_phi)) +
  geom_col() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  # add stars at the right end of each bar
  geom_text(aes(
    label = overdisp_signif,
    y     = overdisp_phi + 0.05 * max_phi
  ),
  hjust = 0) +
  # give a bit of space on the right so stars are visible
  expand_limits(y = max_phi * 1.15) +
  coord_flip() +
  labs(
    title = paste0("Overdispersion of strand proportions: ", plate),
    x     = "Sample",
    y     = "Overdispersion factor (phi)"
  ) +
  theme_bw()

print(p_overdisp)

dev.off()
