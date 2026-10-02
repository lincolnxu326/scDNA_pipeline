#!/usr/bin/env Rscript
# Do any CellenONE droplet measurements separate the read-count gate outcome?
#
# For every well of plate21 and plate22, compare CellenONE-derived metrics across the
# three gate-status groups (PASS / WARN / FAIL). The gate is read-count only, so this
# asks the exploratory question: does anything the camera saw predict whether the cell
# sequenced?
#
# RUN FROM OUTSIDE THE PIPELINE REPO. The repo root carries an .Rprofile that activates
# renv and hijacks .libPaths(), which hides the system ggplot2/ggpubr:
#
#   cd /tmp && Rscript /path/to/scDNA_pipeline/ad_hoc_checks/cellenone_vs_gate.R
#
# Outputs (to --outdir, default ad_hoc_checks/cellenone_vs_gate/):
#   cellenone_vs_gate.pdf    one page per metric, boxplot + jitter, faceted by plate
#   <metric>.png             the same panels individually
#   stats.tsv                every test, effect size and group n
#   wells_merged.tsv         the tidy per-well table the plots were built from

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(ggpubr);  library(stringr)
})

# ---------------------------------------------------------------- configuration ----
DATA_ROOT <- "/nemo/project/proj-tracerX/working/VCAM1_GnT/DATA/384_well"
PLATES    <- c("plate21", "plate22")
SUBPLATES <- 1:4
PASS_CUT  <- 100000   # config.yaml qc_review.usable_reads_pass_cutoff
WARN_CUT  <- 50000    # config.yaml qc_review.usable_reads_warn_cutoff

args   <- commandArgs(trailingOnly = TRUE)
getopt <- function(flag, default) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1]
}
REPO   <- dirname(dirname(normalizePath(sub("--file=", "",
            grep("--file=", commandArgs(FALSE), value = TRUE)[1]))))
OUTDIR <- getopt("--outdir", file.path(REPO, "ad_hoc_checks", "cellenone_vs_gate"))
# Figures go somewhere TRACKED, because docs/CELLENONE_VS_GATE.md embeds them and
# ad_hoc_checks/ is gitignored — a figure the reader cannot see is not a result.
FIGDIR <- getopt("--figdir", file.path(REPO, "docs", "figures", "cellenone_vs_gate"))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
message("Data to:    ", OUTDIR)
message("Figures to: ", FIGDIR)

# ------------------------------------------------------------------ gate status ----
# Recomputed from the BAMs rather than scraped out of review.html, so this analysis
# does not silently inherit whatever the last report happened to be built from.
read_gate <- function(plate) {
  rows <- lapply(SUBPLATES, function(k) {
    sub  <- sprintf("%s_%d", plate, k)
    bdir <- file.path(DATA_ROOT, plate, sub, "bam")
    if (!dir.exists(bdir)) return(NULL)
    wells <- sprintf("W%02d", 1:96)
    reads <- vapply(wells, function(w) {
      f <- file.path(bdir, paste0(w, ".stats.txt"))
      if (!file.exists(f)) return(NA_real_)
      ln <- grep("^SN\treads mapped:", readLines(f, n = 60, warn = FALSE), value = TRUE)
      if (!length(ln)) NA_real_ else as.numeric(strsplit(ln[1], "\t")[[1]][3])
    }, numeric(1))
    tibble(id = paste0(sub, "_", wells), subplate = sub, well = wells, reads = reads)
  })
  bind_rows(rows) %>%
    mutate(status = case_when(is.na(reads)     ~ NA_character_,
                              reads >= PASS_CUT ~ "PASS",
                              reads >= WARN_CUT ~ "WARN",
                              TRUE              ~ "FAIL"))
}

# --------------------------------------------------------------- CellenONE data ----
read_plate <- function(plate) {
  cdir <- file.path(DATA_ROOT, plate, "cellenone")
  wells <- read_tsv(file.path(cdir, "cellenone_wells.tsv"), show_col_types = FALSE,
                    na = c("", "NA"))
  objs  <- read_tsv(file.path(cdir, "objects.tsv"), show_col_types = FALSE,
                    na = c("", "NA"))
  meta  <- jsonlite::fromJSON(file.path(cdir, "run_meta.json"))
  green <- as.numeric(meta$lines$green_x); purple <- as.numeric(meta$lines$purple_x)

  # --- metrics that need deriving from the per-object table ---
  # The isolated cell is the one CellenONE actually dispensed; on both plates it sits
  # left of the green line in 384/384 wells, so this distance is always >= 0.
  iso <- objs %>% filter(category == "isolated") %>%
    group_by(id) %>% slice(1) %>% ungroup() %>%
    transmute(id, dist_iso_to_green = green - x)

  # Nearest object BEYOND the purple line: the one still furthest up the capillary but
  # closest to re-entering the ejection zone. NA when nothing is beyond purple.
  beyond <- objs %>% filter(x > purple) %>%
    group_by(id) %>% summarise(dist_beyond_purple = min(x - purple), .groups = "drop")

  # Any object sitting BETWEEN the lines. Every such object is non-isolated (checked),
  # so this is exactly the contamination-risk indicator.
  between <- objs %>%
    group_by(id) %>%
    summarise(obj_between_lines = any(x > green & x <= purple), .groups = "drop")

  wells %>%
    left_join(iso, by = "id") %>%
    left_join(beyond, by = "id") %>%
    left_join(between, by = "id") %>%
    mutate(plate = plate,
           obj_between_lines = ifelse(is.na(obj_between_lines), FALSE, obj_between_lines),
           green_x = green, purple_x = purple)
}

dat <- lapply(PLATES, function(p) read_plate(p) %>% left_join(read_gate(p), by = "id")) %>%
  bind_rows() %>%
  filter(!is.na(status)) %>%
  mutate(status = factor(status, levels = c("PASS", "WARN", "FAIL")),
         plate  = factor(plate, levels = PLATES))

# ------------------------------------------------------------------ the metrics ----
# label -> column. `rightmost_diameter_um` is the diameter of the object with the
# LARGEST x, i.e. the one furthest from the nozzle and least likely to be ejected —
# it is the object the droplet-call rule keys on, not the isolated cell.
CONT <- c(
  "Isolated cell diameter (um)"            = "diameter_um",
  "Objects detected"                       = "n_objects",
  "Objects in isolation window"            = "n_in_iso_window",
  "Rightmost object diameter (um)"         = "rightmost_diameter_um",
  "Circularity"                            = "circularity",
  "Elongation"                             = "elongation",
  "Blue intensity"                         = "blue_intensity",
  "Orange intensity"                       = "orange_intensity",
  "Isolated cell to green line (px)"       = "dist_iso_to_green",
  "Nearest object past purple line (px)"   = "dist_beyond_purple"
)
BINARY <- c("Object between the lines" = "obj_between_lines")

write_tsv(dat, file.path(OUTDIR, "wells_merged.tsv"))

# ------------------------------------------------------------------------ stats ----
# Kruskal-Wallis across the three groups, then pairwise Wilcoxon (BH). Non-parametric
# throughout: these are counts, bounded ratios and skewed distances, and WARN is tiny.
stat_rows <- list()

for (lab in names(CONT)) {
  col <- CONT[[lab]]
  d <- dat %>% filter(!is.na(.data[[col]])) %>% select(plate, status, value = all_of(col))
  if (nrow(d) < 5 || n_distinct(d$status) < 2) next
  n  <- d %>% count(status)
  kw <- tryCatch(kruskal.test(value ~ status, data = d), error = function(e) NULL)
  pw <- tryCatch(pairwise.wilcox.test(d$value, d$status, p.adjust.method = "BH",
                                      exact = FALSE), error = function(e) NULL)
  getp <- function(a, b) if (is.null(pw)) NA_real_ else
    tryCatch(pw$p.value[a, b], error = function(e) NA_real_)
  # Is pooling safe for THIS metric? If the two plates differ within the FAIL wells
  # (the big, comparable group), then plate is a confounder and the pooled test is
  # mixing two populations rather than gaining power.
  df <- d %>% filter(status == "FAIL")
  pplate <- if (n_distinct(df$plate) == 2)
    tryCatch(wilcox.test(value ~ plate, data = df, exact = FALSE)$p.value,
             error = function(e) NA_real_) else NA_real_
  # Correlation against the raw read count, not just the three-way split. The gate
  # bins a continuous quantity, so a real but modest association could survive
  # binning-induced information loss and show up here instead.
  dc <- dat %>% filter(!is.na(.data[[col]]), !is.na(reads))
  ct <- tryCatch(suppressWarnings(cor.test(dc[[col]], dc$reads, method = "spearman")),
                 error = function(e) NULL)
  med <- d %>% group_by(status) %>% summarise(m = median(value), .groups = "drop")
  gm <- function(s) { v <- med$m[med$status == s]; if (length(v)) v else NA_real_ }
  gn <- function(s) { v <- n$n[n$status == s];    if (length(v)) v else 0L }
  stat_rows[[length(stat_rows) + 1]] <- tibble(
    metric = lab, column = col, test = "Kruskal-Wallis",
    rho_reads = if (is.null(ct)) NA_real_ else unname(ct$estimate),
    p_rho     = if (is.null(ct)) NA_real_ else ct$p.value,
    p_omnibus = if (is.null(kw)) NA_real_ else kw$p.value,
    p_PASS_vs_FAIL = getp("FAIL", "PASS"),
    p_PASS_vs_WARN = getp("WARN", "PASS"),
    p_WARN_vs_FAIL = getp("WARN", "FAIL"),
    p_plate_effect = pplate,
    n_PASS = gn("PASS"), n_WARN = gn("WARN"), n_FAIL = gn("FAIL"),
    median_PASS = gm("PASS"), median_WARN = gm("WARN"), median_FAIL = gm("FAIL"))
}

for (lab in names(BINARY)) {
  col <- BINARY[[lab]]
  d  <- dat %>% select(plate, status, value = all_of(col))
  tb <- table(d$status, d$value)
  if (nrow(tb) < 2 || ncol(tb) < 2) next
  ft <- tryCatch(fisher.test(tb), error = function(e) NULL)
  pplate <- tryCatch(fisher.test(table(d$plate[d$status == "FAIL"],
                                       d$value[d$status == "FAIL"]))$p.value,
                     error = function(e) NA_real_)
  prop <- d %>% group_by(status) %>% summarise(p = mean(value), n = n(), .groups = "drop")
  gp <- function(s) { v <- prop$p[prop$status == s]; if (length(v)) v else NA_real_ }
  gn <- function(s) { v <- prop$n[prop$status == s]; if (length(v)) v else 0L }
  dcb <- dat %>% filter(!is.na(reads))
  ctb <- tryCatch(suppressWarnings(cor.test(as.numeric(dcb[[col]]), dcb$reads,
                                            method = "spearman")),
                  error = function(e) NULL)
  stat_rows[[length(stat_rows) + 1]] <- tibble(
    metric = lab, column = col, test = "Fisher exact",
    rho_reads = if (is.null(ctb)) NA_real_ else unname(ctb$estimate),
    p_rho     = if (is.null(ctb)) NA_real_ else ctb$p.value,
    p_omnibus = if (is.null(ft)) NA_real_ else ft$p.value,
    p_PASS_vs_FAIL = NA_real_, p_PASS_vs_WARN = NA_real_, p_WARN_vs_FAIL = NA_real_,
    p_plate_effect = pplate,
    n_PASS = gn("PASS"), n_WARN = gn("WARN"), n_FAIL = gn("FAIL"),
    median_PASS = gp("PASS"), median_WARN = gp("WARN"), median_FAIL = gp("FAIL"))
}

stats <- bind_rows(stat_rows) %>%
  mutate(p_omnibus_BH = p.adjust(p_omnibus, "BH"),
         p_rho_BH      = p.adjust(p_rho, "BH")) %>%
  arrange(p_rho)
write_tsv(stats, file.path(OUTDIR, "stats.tsv"))

# ------------------------------------------------------------------------ plots ----
FILL  <- c(PASS = "#2e7d4f", WARN = "#2f6cad", FAIL = "#b23a2f")
PSHAPE <- c(plate21 = 16, plate22 = 17)

panel <- function(lab, col) {
  d  <- dat %>% filter(!is.na(.data[[col]])) %>% select(plate, status, value = all_of(col))
  ns <- d %>% count(status) %>% mutate(lbl = paste0("n=", n))
  ggplot(d, aes(status, value)) +
    geom_boxplot(aes(fill = status), outlier.shape = NA, width = .5, alpha = .35,
                 colour = "grey25") +
    geom_jitter(aes(colour = status, shape = plate), width = .2, height = 0,
                size = 1.2, alpha = .5) +
    stat_compare_means(method = "kruskal.test", label = "p.format",
                       size = 3.4, label.y.npc = .97) +
    geom_text(data = ns, aes(x = status, y = -Inf, label = lbl), inherit.aes = FALSE,
              vjust = -0.6, size = 3, colour = "grey30") +
    scale_fill_manual(values = FILL) + scale_colour_manual(values = FILL) +
    scale_shape_manual(values = PSHAPE, name = NULL) +
    guides(fill = "none", colour = "none") +
    labs(title = lab, subtitle = "plate21 + plate22 pooled", x = NULL, y = lab) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          legend.position = "top", legend.justification = "right")
}

binary_panel <- function(lab, col) {
  d <- dat %>% select(status, value = all_of(col)) %>%
    group_by(status) %>% summarise(prop = mean(value), n = n(), .groups = "drop")
  pv <- stats$p_omnibus[stats$metric == lab][1]
  ggplot(d, aes(status, prop, fill = status)) +
    geom_col(width = .55, alpha = .8) +
    geom_text(aes(label = paste0(round(prop * 100, 1), "%\nn=", n)),
              vjust = -0.25, size = 3, colour = "grey20") +
    annotate("text", x = 2, y = Inf, vjust = 1.6, size = 3.4,
             label = paste0("Fisher p = ", signif(pv, 3))) +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .25))) +
    scale_fill_manual(values = FILL) +
    labs(title = lab, subtitle = "plate21 + plate22 pooled", x = NULL,
         y = "wells with such an object") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40", size = 9))
}

scatter <- function(lab, col) {
  d <- dat %>% filter(!is.na(.data[[col]]), !is.na(reads), reads > 0) %>%
    mutate(value = as.numeric(.data[[col]]))
  r <- stats %>% filter(metric == lab) %>% slice(1)
  ggplot(d, aes(value, reads)) +
    geom_point(aes(colour = status, shape = plate), size = 1.2, alpha = .55) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
                colour = "grey25", linewidth = .6) +
    geom_hline(yintercept = c(WARN_CUT, PASS_CUT), linetype = "dashed",
               colour = "grey55", linewidth = .35) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.6, size = 3.4,
             label = sprintf("Spearman rho = %.3f, p = %.3g", r$rho_reads, r$p_rho)) +
    scale_y_log10(labels = scales::comma) +
    scale_colour_manual(values = FILL) + scale_shape_manual(values = PSHAPE, name = NULL) +
    guides(colour = "none") +
    labs(title = paste0(lab, " vs read count"),
         subtitle = "plate21 + plate22 pooled; dashed lines = the 50k / 100k gate",
         x = lab, y = "usable reads (log10)") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey40", size = 9),
          legend.position = "top", legend.justification = "right")
}

plots <- c(lapply(names(CONT),   function(l) panel(l, CONT[[l]])),
           lapply(names(BINARY), function(l) binary_panel(l, BINARY[[l]])),
           lapply(names(CONT),   function(l) scatter(l, CONT[[l]])))
names(plots) <- c(names(CONT), names(BINARY), paste0(names(CONT), " vs reads"))

pdf(file.path(OUTDIR, "cellenone_vs_gate.pdf"), width = 7.5, height = 4.6, onefile = TRUE)
for (p in plots) print(p)
invisible(dev.off())

for (nm in names(plots)) {
  f <- file.path(FIGDIR, paste0(str_replace_all(str_to_lower(nm), "[^a-z0-9]+", "_"), ".png"))
  ggsave(f, plots[[nm]], width = 7.5, height = 4.6, dpi = 150)
}
# the summary table travels with the figures so the doc's numbers stay auditable
write_tsv(stats, file.path(FIGDIR, "stats.tsv"))

# ---------------------------------------------------------------------- summary ----
message("\n", strrep("-", 78))
message("Kruskal-Wallis / Fisher across PASS / WARN / FAIL, most significant first")
message(strrep("-", 78))
stats %>%
  transmute(metric,
            p = signif(p_omnibus, 3), p_BH = signif(p_omnibus_BH, 3),
            `PASS vs FAIL` = signif(p_PASS_vs_FAIL, 3),
            rho_reads = signif(rho_reads, 3), `p(rho)` = signif(p_rho, 3),
            `plate confound` = signif(p_plate_effect, 3),
            n = paste(n_PASS, n_WARN, n_FAIL, sep = "/"),
            medians = paste(signif(median_PASS, 3), signif(median_WARN, 3),
                            signif(median_FAIL, 3), sep = " / ")) %>%
  as.data.frame() %>% print(row.names = FALSE)
message("\nn and medians are PASS / WARN / FAIL.")
message("Both plates pooled. WARN is n=9 in total — underpowered; read PASS vs FAIL.")
message("`plate confound` = plate21 vs plate22 WITHIN the FAIL wells. Small values mean ",
        "the plates differ for that metric, so the pooled test mixes two populations.")
message("Done: ", OUTDIR)
